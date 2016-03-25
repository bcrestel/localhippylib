# Copyright (c) 2016, The University of Texas at Austin & University of
# California, Merced.
#
# All Rights reserved.
# See file COPYRIGHT for details.
#
# This file is part of the hIPPYlib library. For more information and source code
# availability see https://hippylib.github.io.
#
# hIPPYlib is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License (as published by the Free
# Software Foundation) version 2.1 dated February 1999.

from dolfin import *
import numpy as np
from timeDependentVector import *
from cgsolverSteihaug import *

STATE = 0
PARAMETER = 1
ADJOINT = 2

class ReducedHessian:
    """
    This class implements matrix free application of the reduced hessian operator.
    The constructor takes the following parameters:
    - model:               the object which contains the description of the problem.
    - innerTol:            the relative tolerance for the solution of the incremental
                           forward and adjoint problems.
    """
    def __init__(self, model, innerTol):
        """
        Construct the reduced Hessian Operator:
        """
        self.model = model
        self.tol = innerTol
        self.ncalls = 0
        
        self.rhs_fwd = model.generate_vector(STATE)
        self.rhs_adj = model.generate_vector(ADJOINT)
        self.rhs_adj2 = model.generate_vector(ADJOINT)
        self.uhat    = model.generate_vector(STATE)
        self.phat    = model.generate_vector(ADJOINT)
        self.yhelp = model.generate_vector(PARAMETER)
    
    def init_vector(self, x, dim):
        """
        Reshape the Vector x so that it is compatible with the reduced Hessian
        operator.
        Parameters:
        - x: the vector to reshape
        - dim: if 0 then x will be reshaped to be compatible with the range of
               the reduced Hessian
               if 1 then x will be reshaped to be compatible with the domain of
               the reduced Hessian
               
        Note: Since the reduced Hessian is a self adjoint operator, the range and
              the domain is the same. Either way, we choosed to add the parameter
              dim for consistency with the interface of Matrix in dolfin.
        """
        self.model.init_parameter(x)
        
    def mult(self,x,y):
        """
        Apply the Gauss Newton approximation of the reduced Hessian to the vector x
        Return the result in y.        
        """
        self.model.applyC(x, self.rhs_fwd)
        self.model.solveFwdIncremental(self.uhat, self.rhs_fwd, self.tol)
        self.model.applyWuu(self.uhat, self.rhs_adj)
        self.model.solveAdjIncremental(self.phat, self.rhs_adj, self.tol)
        self.model.applyCt(self.phat, y)
        
        self.model.applyR(x,self.yhelp)
        y.axpy(1., self.yhelp)


class LaplacianRegularization:
    def __init__(self, Vh, gamma, delta, a0=None, rel_tol=1e-12, max_iter=100):
        """
        Construct the Regularization.
        Input:
        - Vh:              the finite element space for the parameter
        - gamma and delta: the coefficient in the PDE
        - Theta:           the s.p.d. tensor for anisotropic diffusion of the pde
        - a0:              the reference values
        """        
        assert delta != 0., "Singular regularization are not supported"
        self.Vh = Vh
        
        trial = TrialFunction(Vh)
        test  = TestFunction(Vh)
        
        varfL = gamma*inner(nabla_grad(trial), nabla_grad(test))*dx
        varfM = delta*inner(trial,test)*dx
        
        self.M = assemble(varfM)
        self.R = assemble(varfL + varfM)
        
        self.Rsolver = PETScKrylovSolver("cg", "amg")
        self.Rsolver.set_operator(self.R)
        self.Rsolver.parameters["maximum_iterations"] = max_iter
        self.Rsolver.parameters["relative_tolerance"] = rel_tol
        self.Rsolver.parameters["error_on_nonconvergence"] = True
        self.Rsolver.parameters["nonzero_initial_guess"] = False

        self.a0 = a0
        
        if self.a0 is None:
            self.a0 = Vector()
            self.init_vector(self.a0, 0)
        
    def init_vector(self,x,dim):
        """
        Inizialize a vector x to be compatible with the range/domain of R.
        If dim == "noise" inizialize x to be compatible with the size of
        white noise used for sampling.
        """
        self.R.init_vector(x,dim)
            


class TimeDependentAD:    
    def __init__(self, mesh, Vh, t_init, t_final, t_1, dt, wind_velocity, gls_stab, reg):
        self.mesh = mesh
        self.Vh = Vh
        self.t_init = t_init
        self.t_final = t_final
        self.t_1 = t_1
        self.dt = dt
        self.sim_times = np.arange(self.t_init, self.t_final+.5*self.dt, self.dt)
        
        u = TrialFunction(Vh[STATE])
        v = TestFunction(Vh[STATE])
        
        kappa = Constant(.001)
        dt_expr = Constant(self.dt)
        
        r_trial = u + dt_expr*( -div(kappa*nabla_grad(u))+ inner(wind_velocity, nabla_grad(u)) )
        r_test  = v + dt_expr*( -div(kappa*nabla_grad(v))+ inner(wind_velocity, nabla_grad(v)) )

        
        h = CellSize(mesh)
        vnorm = sqrt(inner(wind_velocity, wind_velocity))
        if gls_stab:
            tau = Min((h*h)/(Constant(2.)*kappa), h/vnorm )
        else:
            tau = Constant(0.)
                            
        self.M = assemble( inner(u,v)*dx )
        self.M_stab = assemble( inner(u, v+tau*r_test)*dx )
        self.Mt_stab = assemble( inner(u+tau*r_trial,v)*dx )
        Nvarf  = (inner(kappa *nabla_grad(u), nabla_grad(v)) + inner(wind_velocity, nabla_grad(u))*v )*dx
        Ntvarf  = (inner(kappa *nabla_grad(v), nabla_grad(u)) + inner(wind_velocity, nabla_grad(v))*u )*dx
        self.N  = assemble( Nvarf )
        self.Nt = assemble(Ntvarf)
        stab = assemble( tau*inner(r_trial, r_test)*dx)
        self.L = self.M + dt*self.N + stab
        self.Lt = self.M + dt*self.Nt + stab
        
        boundaries = FacetFunction("size_t", mesh)
        boundaries.set_all(0)

        class InsideBoundary(SubDomain):
            def inside(self,x,on_boundary):
                x_in = x[0] > DOLFIN_EPS and x[0] < 1 - DOLFIN_EPS
                y_in = x[1] > DOLFIN_EPS and x[1] < 1 - DOLFIN_EPS
                return on_boundary and x_in and y_in
            
        Gamma_M = InsideBoundary()
        Gamma_M.mark(boundaries,1)
        ds_marked = Measure("ds")[boundaries]
        
        self.Q = assemble( self.dt*inner(u, v) * ds_marked(1) )

        self.reg = reg
        
        self.solver = PETScKrylovSolver("gmres", "ilu")
        self.solver.set_operator(self.L)
        
        self.solvert = PETScKrylovSolver("gmres", "ilu")
        self.solvert.set_operator(self.Lt)
                        
        self.ud = self.generate_vector(STATE)
        self.noise_variance = 0
                
    def generate_vector(self, component = "ALL"):
        if component == "ALL":
            u = TimeDependentVector(self.sim_times)
            u.initialize(self.Q, 0)
            a = Vector()
            self.reg.init_vector(a,0)
            p = TimeDependentVector(self.sim_times)
            p.initialize(self.Q, 0)
            return [u, a, p]
        elif component == STATE:
            u = TimeDependentVector(self.sim_times)
            u.initialize(self.Q, 0)
            return u
        elif component == PARAMETER:
            a = Vector()
            self.reg.init_vector(a,0)
            return a
        elif component == ADJOINT:
            p = TimeDependentVector(self.sim_times)
            p.initialize(self.Q, 0)
            return p
        else:
            raise
    
    def init_parameter(self, a):
        self.reg.init_vector(a,0)
          
    def cost(self, x):
        Rdx = Vector()
        self.reg.init_vector(Rdx,0)
        dx = x[PARAMETER] - self.reg.a0
        self.reg.R.mult(dx, Rdx)
        reg = .5*Rdx.inner(dx)
        
        u  = Vector()
        ud = Vector()
        self.Q.init_vector(u,0)
        self.Q.init_vector(ud,0)
    
        misfit = 0
        for t in np.arange(self.t_1, self.t_final+(.5*self.dt), self.dt):
            x[STATE].retrieve(u,t)
            self.ud.retrieve(ud,t)
            diff = u - ud
            Qdiff = self.Q * diff
            misfit += .5/self.noise_variance*Qdiff.inner(diff)
            
        c = misfit + reg
                
        return [c, reg, misfit]
    
    def solveFwd(self, out, x, tol=1e-9):
        out.zero()
        uold = x[PARAMETER]
        u = Vector()
        rhs = Vector()
        self.M.init_vector(rhs, 0)
        self.M.init_vector(u, 0)
        self.solver.parameters["relative_tolerance"] = tol
        t = self.t_init
        while t < self.t_final:
            t += self.dt
            self.M_stab.mult(uold, rhs)
            self.solver.solve(u, rhs)
            out.store(u,t)
            uold = u
    
    def solveAdj(self, out, x, tol=1e-9):
        out.zero()
        pold = Vector()
        self.M.init_vector(pold,0)    
        p = Vector()
        self.M.init_vector(p,0)
        rhs = Vector()
        self.M.init_vector(rhs,0) 
        rhs_obs = Vector()
        
        u = Vector()
        self.M.init_vector(u,0)
        ud = Vector()
        self.M.init_vector(ud,0)
  
        self.solvert.parameters["relative_tolerance"] = tol
                
        t = self.t_final
        while t > self.t_init:
            self.Mt_stab.mult(pold,rhs)
            if t > self.t_1 - .5*self.dt:
                x[STATE].retrieve(u,t)
                self.ud.retrieve(ud,t)
                ud.axpy(-1., u)
                self.Q.mult(ud,rhs_obs)
#                print "t = ", t, "solveAdj ||ud-u||_inf = ", ud.norm("linf"), " ||rhs_obs|| = ", rhs_obs.norm("linf")
                rhs.axpy(1./self.noise_variance, rhs_obs)
                
            self.solvert.solve(p, rhs)
            pold = p
            out.store(p, t)
            t -= self.dt
            
            
            
    def evalGradientParameter(self,x, mg):
        self.reg.init_vector(mg,1)
        dx = x[PARAMETER] - self.reg.a0
        self.reg.R.mult(dx, mg)
        
        p0 = Vector()
        self.Q.init_vector(p0,0)
        x[ADJOINT].retrieve(p0, self.t_init + self.dt)
        
        mg.axpy(-1., self.Mt_stab*p0)
        
        g = Vector()
        self.M.init_vector(g,1)
        
        s = PETScKrylovSolver("cg", "jacobi")
        s.parameters["relative_tolerance"] = 1e-9
        s.set_operator(self.M)
        s.solve(g,mg)
        
        
        grad_norm = g.inner(mg)
        
        return grad_norm
        
            
    def solveFwdIncremental(self, sol, rhs, tol):
        sol.zero()
        uold = Vector()
        u = Vector()
        Muold = Vector()
        myrhs = Vector()
        self.M.init_vector(uold, 0)
        self.M.init_vector(u, 0)
        self.M.init_vector(Muold, 0)
        self.M.init_vector(myrhs, 0)
        self.solver.parameters["relative_tolerance"] = tol
        t = self.t_init
        while t < self.t_final:
            t += self.dt
            self.M_stab.mult(uold, Muold)
            rhs.retrieve(myrhs, t)
            myrhs.axpy(1., Muold)
            self.solver.solve(u, myrhs)
            sol.store(u,t)
            uold = u


        
    def solveAdjIncremental(self, sol, rhs, tol):
        sol.zero()
        pold = Vector()
        p = Vector()
        Mpold = Vector()
        myrhs = Vector()
        self.M.init_vector(pold, 0)
        self.M.init_vector(p, 0)
        self.M.init_vector(Mpold, 0)
        self.M.init_vector(myrhs, 0)
        self.solvert.parameters["relative_tolerance"] = tol
        t = self.t_final
        while t > self.t_init:
            self.Mt_stab.mult(pold, Mpold)
            rhs.retrieve(myrhs, t)
            myrhs.axpy(1., Mpold)
            self.solvert.solve(p, myrhs)
            sol.store(p,t)
            pold = p
            t -= self.dt
    
    def applyC(self, da, out):
        out.zero()
        myout = Vector()
        self.M.init_vector(myout, 0)
        self.M_stab.mult(da,myout)
        myout *= -1.
        t = self.t_init + self.dt
        out.store(myout,t)
        
        myout.zero()
        while t < self.t_final:
            t += self.dt
            out.store(myout,t)
    
    def applyCt(self, dp, out):
        t = self.t_init + self.dt
        dp0 = Vector()
        self.M.init_vector(dp0,0)
        dp.retrieve(dp0, t)
        dp0 *= -1.
        self.Mt_stab.mult(dp0, out)

    
    def applyWuu(self, du, out):
        out.zero()
        myout = Vector()
        self.Q.init_vector(myout,0)
        myout.zero()
        
        t = self.t_init + self.dt
        while t < self.t_1 - .5*self.dt:
            out.store(myout, t)
            t += self.dt
            
        mydu  = Vector()
        self.Q.init_vector(mydu,0)
        while t < self.t_final+(.5*self.dt):
            du.retrieve(mydu,t)
            self.Q.mult(mydu, myout)
            myout *= 1./self.noise_variance
            out.store(myout, t)
            t += self.dt
        
    def applyR(self, da, out):
        self.reg.R.mult(da,out)
          
    def exportState(self, x, filename, varname):
        out_file = File(filename)
        ufunc = Function(self.Vh[STATE], name=varname)
        t = self.t_init
        out_file << (Function(self.Vh[STATE], x[PARAMETER], name=varname),t)
        while t < self.t_final:
            t += self.dt
            x[STATE].retrieve(ufunc.vector(), t)
            out_file << (ufunc, t)
            
        
        
def v_boundary(x,on_boundary):
    return on_boundary

def q_boundary(x,on_boundary):
    return x[0] < DOLFIN_EPS and x[1] < DOLFIN_EPS
        
def computeVelocityField(mesh):
    Xh = VectorFunctionSpace(mesh,'Lagrange', 2)
    Wh = FunctionSpace(mesh, 'Lagrange', 1)
    XW = MixedFunctionSpace([Xh, Wh])

    
    Re = 1e2
    
    g = Expression(('0.0','(x[0] < 1e-14) - (x[0] > 1 - 1e-14)'))
    bc1 = DirichletBC(XW.sub(0), g, v_boundary)
    bc2 = DirichletBC(XW.sub(1), Constant(0), q_boundary, 'pointwise')
    bcs = [bc1, bc2]
    
    vq = Function(XW)
    (v,q) = split(vq)
    (v_test, q_test) = TestFunctions (XW)
    
    def strain(v):
        return sym(nabla_grad(v))
    
    F = ( (2./Re)*inner(strain(v),strain(v_test))+ inner (nabla_grad(v)*v, v_test)
           - (q * div(v_test)) + ( div(v) * q_test) ) * dx
           
    solve(F == 0, vq, bcs, solver_parameters={"newton_solver":
                                         {"relative_tolerance":1e-4, "maximum_iterations":100}})
        
    return v
    

        
if __name__ == "__main__":
    set_log_active(False)
    np.random.seed(1)
    sep = "\n"+"#"*80+"\n"
    print sep, "Set up the mesh and finite element spaces.\n","Compute wind velocity", sep
    mesh = refine( Mesh("ad_20.xml") )
    wind_velocity = computeVelocityField(mesh)
    Vh = FunctionSpace(mesh, "Lagrange", 2)
    print "Number of dofs: {0}".format( Vh.dim() )
    
    print sep, "Set up Prior Information and model", sep
    
    true_initial_condition = interpolate(Expression('min(0.5,exp(-100*(pow(x[0]-0.35,2) +  pow(x[1]-0.7,2))))'), Vh).vector()

    gamma = 1
    delta = 1e1
    reg = LaplacianRegularization(Vh, gamma, delta)
    
    reg.a0 = interpolate(Expression('0.5'), Vh).vector()
    
    print "Regularization: (delta - gamma*Laplacian): delta={0}, gamma={1}".format(delta, gamma)

    problem = TimeDependentAD(mesh, [Vh,Vh,Vh], 0., 4., 1., .2, wind_velocity, True, reg)
    
    print sep, "Generate synthetic observation", sep
    rel_noise = 0.001
    utrue = problem.generate_vector(STATE)
    x = [utrue, true_initial_condition, None]
    problem.solveFwd(x[STATE], x, 1e-9)
    MAX = utrue.norm("linf", "linf")
    noise_std_dev = rel_noise * MAX
    problem.ud.copy(utrue)
    problem.ud.randn_perturb(noise_std_dev)
    problem.noise_variance = noise_std_dev*noise_std_dev
    
    print sep, "Compute the reduced gradient and hessian", sep
    [u,a,p] = problem.generate_vector()
    problem.solveFwd(u, [u,a,p], 1e-12)
    problem.solveAdj(p, [u,a,p], 1e-12)
    mg = problem.generate_vector(PARAMETER)
    grad_norm = problem.evalGradientParameter([u,a,p], mg)
        
    print "(g,g) = ", grad_norm
    
    H = ReducedHessian(problem, 1e-12) 
    
    print sep, "Find the MAP point", sep
            
    solver = CGSolverSteihaug()
    solver.set_operator(H)
    solver.set_preconditioner( reg.Rsolver )
    solver.parameters["print_level"] = 1
    solver.parameters["rel_tolerance"] = 1e-6
    solver.solve(a, -mg)
    problem.solveFwd(u, [u,a,p], 1e-12)
 
    total_cost, reg_cost, misfit_cost = problem.cost([u,a,p])
    print "Total cost {0:5g}; Reg Cost {1:5g}; Misfit {2:5g}".format(total_cost, reg_cost, misfit_cost)
        
    print sep, "Save results", sep  
    problem.exportState([u,a,p], "results/conc.pvd", "concentration")
    problem.exportState([utrue,true_initial_condition,p], "results/true_conc.pvd", "concentration")
    problem.exportState([problem.ud,true_initial_condition,p], "results/noisy_conc.pvd", "concentration")

    print sep, "Visualize results", sep 
    plot(Function(Vh,a, name = "Initial Condition"))
    
