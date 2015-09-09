from dolfin import *
import numpy as np
import sys
sys.path.append( "../" )
from pylib import *
import matplotlib.pyplot as plt

class TimeDependentAD:    
    def __init__(self, mesh, Vh, t_init, t_final, t_1, dt, wind_velocity, true_initial_condition, gls_stab, Prior):
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
        
        boundary_parts = MeshFunction("size_t",mesh, mesh.topology().dim()-1)

        class InsideBoundary(SubDomain):
            def inside(self,x,on_boundary):
                x_in = x[0] > DOLFIN_EPS and x[0] < 1 - DOLFIN_EPS
                y_in = x[1] > DOLFIN_EPS and x[1] < 1 - DOLFIN_EPS
                return on_boundary and x_in and y_in
            
        Gamma_M = InsideBoundary()
        Gamma_M.mark(boundary_parts,1)
        
        self.Q = assemble( self.dt*inner(u, v) * ds(1), exterior_facet_domains=boundary_parts )

        self.Prior = Prior
        
        self.solver = PETScKrylovSolver("gmres", "ilu")
        self.solver.set_operator(self.L)
        
        self.solvert = PETScKrylovSolver("gmres", "ilu")
        self.solvert.set_operator(self.Lt)
                        
        self.true_init = true_initial_condition.copy()

        self.ud = self.generate_vector(STATE)
        self.noise_var = self.computeObservation(self.ud)
        print "Variance of Noise: ", self.noise_var
                
    def generate_vector(self, component = "ALL"):
        if component == "ALL":
            u = TimeDependentVector(self.sim_times)
            u.initialize(self.Q, 0)
            control = Vector()
            self.Prior.init_vector(control,0)
            p = TimeDependentVector(self.sim_times)
            p.initialize(self.Q, 0)
            return [u, control, p]
        elif component == STATE:
            u = TimeDependentVector(self.sim_times)
            u.initialize(self.Q, 0)
            return u
        elif component == CONTROL:
            control = Vector()
            self.Prior.init_vector(control,0)
            return control
        elif component == ADJOINT:
            p = TimeDependentVector(self.sim_times)
            p.initialize(self.Q, 0)
            return p
        else:
            raise
    
    def init_control(self, a):
        self.Prior.init_vector(a,0)
        
    def getIdentityMatrix(self, component):
        Xh = self.Vh[component]
        test = TestFunction(Xh)
        trial = TrialFunction(Xh)
        
        I = assemble(test*trial*dx)
        I.zero()
        I.ident_zeros()
        
        return I
        
          
    def cost(self, x):
        Rdx = Vector()
        self.Prior.init_vector(Rdx,0)
        dx = x[CONTROL] - self.Prior.mean
        self.Prior.R.mult(dx, Rdx)
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
            misfit += .5/self.noise_var*Qdiff.inner(diff)
            
        c = misfit + reg
                
        return [c, reg, misfit]
    
    def solveFwd(self, out, x, tol=1e-9):
        out.zero()
        uold = x[CONTROL]
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

    def computeObservation(self, ud):
        ud.zero()
        uold_func = Function(self.Vh[CONTROL], self.true_init)
        uold = uold_func.vector()
        
        np.random.seed(1)
         
        u = Vector()
        rhs = Vector()
        self.M.init_vector(rhs, 0)
        self.M.init_vector(u, 0)
        self.solver.parameters["relative_tolerance"] = 1e-9
        t = self.t_init
        MAX = uold.norm("linf")
        while t < self.t_final:
            t += self.dt
            self.M_stab.mult(uold, rhs)
            self.solver.solve(u, rhs)
            uold.set_local(u.array())
            umax = u.norm("linf")
            if umax > MAX:
                MAX = umax
            ud.store(u,t)
        
        std_dev = 0.005*MAX
        
        out_file = File("data/conc.pvd")
        while t < self.t_final:
            t += self.dt
            self.ud.retrieve(u,t)
            noise = std_dev * np.random.normal(0, 1, len(u.array()))
            u.set_local(u.array() + noise)
            out_file << (Function(self.Vh[STATE],u, name='data_conc'),t)
            ud.store(u,t)
            
        return std_dev*std_dev

    
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
                rhs.axpy(1./self.noise_var, rhs_obs)
                
            self.solvert.solve(p, rhs)
            pold = p
            out.store(p, t)
            t -= self.dt
            
            
            
    def evalGradientControl(self,x, mg):
        self.Prior.init_vector(mg,1)
        dx = x[CONTROL] - self.Prior.mean
        self.Prior.R.mult(dx, mg)
        
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
        
    
    def setPointForHessianEvaluations(self, x):
        """
        Specify the point x = [u,a,p] at which the Hessian operator (or the Gauss-Newton approximation)
        need to be evaluated.
        
        Nothing to do since the problem is linear
        """
        return

        
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
            myout *= 1./self.noise_var
            out.store(myout, t)
            t += self.dt
    
    def applyWua(self, da, out):
        out.zero()

    
    def applyWau(self, du, out):
        out.zero()
    
    def applyR(self, da, out):
        self.Prior.R.mult(da,out)
    
    def applyRaa(self, da, out):
        out.zero()
        
    def exportState(self, x, filename, varname):
        out_file = File(filename)
        ufunc = Function(self.Vh[STATE], name=varname)
        t = self.t_init
        out_file << (Function(self.Vh[STATE], x[CONTROL], name=varname),t)
        while t < self.t_final:
            t += self.dt
            x[STATE].retrieve(ufunc.vector(), t)
            out_file << (ufunc, t)

