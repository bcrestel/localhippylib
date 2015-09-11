from dolfin import *
import sys
sys.path.append( "../../" )
from hippylib import *
import numpy as np
import matplotlib.pyplot as plt

def u_boundary(x, on_boundary):
    return on_boundary

class Poisson:
    def __init__(self, mesh, Vh, Prior):
        """
        Construct a model by proving
        - the mesh
        - the finite element spaces for the STATE/ADJOINT variable and the control variable
        - the Prior information
        """
        self.mesh = mesh
        self.Vh = Vh
        
        # Initialize Expressions
        self.atrue = Expression('log(2 + 7*(pow(pow(x[0] - 0.5,2) + pow(x[1] - 0.5,2),0.5) > 0.2))')
        self.f = Expression("1.0")
        self.u_o = Vector()
        
        self.u_bdr = Expression("0.0")
        self.u_bdr0 = Expression("0.0")
        self.bc = DirichletBC(self.Vh[STATE], self.u_bdr, u_boundary)
        self.bc0 = DirichletBC(self.Vh[STATE], self.u_bdr0, u_boundary)
        
        # Assemble constant matrices      
        self.M = self.assembleM()
        self.Prior = Prior
        self.Wuu = self.assembleWuu()
        

        self.computeObservation(self.u_o)
                
        self.A = []
        self.At = []
        self.C = []
        self.Raa = []
        self.Wau = []
        
    def generate_vector(self, component="ALL"):
        """
        Return the list x=[u,a,p] where:
        - u is any object that describes the state variable
        - a is a Vector object that describes the control variable.
          (Need to support linear algebra operations)
        - p is any object that describes the adjoint variable
        
        If component is STATE, CONTROL, or ADJOINT return x[component]
        """
        if component == "ALL":
            x = [Vector(), Vector(), Vector()]
            self.Wuu.init_vector(x[STATE],0)
            self.Prior.init_vector(x[CONTROL],0)
            self.Wuu.init_vector(x[ADJOINT], 0)
        elif component == STATE:
            x = Vector()
            self.Wuu.init_vector(x,0)
        elif component == CONTROL:
            x = Vector()
            self.Prior.init_vector(x,0)
        elif component == ADJOINT:
            x = Vector()
            self.Wuu.init_vector(x,0)
            
        return x
    
    def init_control(self, a):
        """
        Reshape a so that it is compatible with the control variable
        """
        self.Prior.init_vector(a,0)
        
    def assembleA(self,x, assemble_adjoint = False, assemble_rhs = False):
        """
        Assemble the matrices and rhs for the forward/adjoint problems
        """
        trial = TrialFunction(self.Vh[STATE])
        test = TestFunction(self.Vh[STATE])
        c = Function(self.Vh[CONTROL], x[CONTROL])
        Avarf = inner(exp(c)*nabla_grad(trial), nabla_grad(test))*dx
        if not assemble_adjoint:
            bform = inner(self.f, test)*dx
            Matrix, rhs = assemble_system(Avarf, bform, self.bc)
        else:
            # Assemble the adjoint of A (i.e. the transpose of A)
            s = Function(self.Vh[STATE], x[STATE])
            obs = Function(self.Vh[STATE], self.u_o)
            bform = inner(obs - s, test)*dx
            Matrix, rhs = assemble_system(adjoint(Avarf), bform, self.bc0)
            
        if assemble_rhs:
            return Matrix, rhs
        else:
            return Matrix
    
    def assembleC(self, x):
        """
        Assemble the derivative of the forward problem with respect to the control
        """
        trial = TrialFunction(self.Vh[CONTROL])
        test = TestFunction(self.Vh[STATE])
        s = Function(Vh[STATE], x[STATE])
        c = Function(Vh[CONTROL], x[CONTROL])
        Cvarf = inner(exp(c) * trial * nabla_grad(s), nabla_grad(test)) * dx
        C = assemble(Cvarf)
#        print "||c||", x[CONTROL].norm("l2"), "||s||", x[STATE].norm("l2"), "||C||", C.norm("linf")
        self.bc0.zero(C)
        return C
   
    def assembleM(self):
        """
        Assemble the mass matrix in the control space.
        This is needed in evalGradientControl to compute the L2 norm of the gradient
        """
        trial = TrialFunction(self.Vh[CONTROL])
        test  = TestFunction(self.Vh[CONTROL])
        varf = inner(trial, test)*dx
        return assemble(varf)
    
    def assembleWuu(self):
        """
        Assemble the misfit operator
        """
        trial = TrialFunction(self.Vh[STATE])
        test = TestFunction(self.Vh[STATE])
        varf = inner(trial, test)*dx
        Wuu = assemble(varf)
        dummy = Vector()
        Wuu.init_vector(dummy,0)
        self.bc0.zero_columns(Wuu, dummy)
        self.bc0.zero(Wuu)
        return Wuu
    
    def assembleWau(self, x):
        """
        Assemble the derivative of the control equation with respect to the state
        """
        trial = TrialFunction(self.Vh[STATE])
        test  = TestFunction(self.Vh[CONTROL])
        a = Function(self.Vh[ADJOINT], x[ADJOINT])
        c = Function(self.Vh[CONTROL], x[CONTROL])
        varf = inner(exp(c)*nabla_grad(trial),nabla_grad(a))*test*dx
        Wau = assemble(varf)
        dummy = Vector()
        Wau.init_vector(dummy,0)
        self.bc0.zero_columns(Wau, dummy)
        return Wau
    
    def assembleRaa(self, x):
        """
        Assemble the derivative of the control equation with respect to the control (Newton method)
        """
        trial = TrialFunction(self.Vh[CONTROL])
        test  = TestFunction(self.Vh[CONTROL])
        s = Function(self.Vh[STATE], x[STATE])
        c = Function(self.Vh[CONTROL], x[CONTROL])
        a = Function(self.Vh[ADJOINT], x[ADJOINT])
        varf = inner(nabla_grad(a),exp(c)*nabla_grad(s))*trial*test*dx
        return assemble(varf)

        
    def computeObservation(self, u_o):
        """
        Compute the syntetic observation
        """
        at = interpolate(self.atrue, Vh[CONTROL])
        x = [self.generate_vector(STATE), at.vector(), None]
        A, b = self.assembleA(x, assemble_rhs = True)
        
        A.init_vector(u_o, 1)
        solve(A, u_o, b)
        
        # Create noisy data, ud
        MAX = u_o.norm("linf")
        noise = .01 * MAX * np.random.normal(0, 1, len(u_o.array()))
        u_o.set_local(u_o.array() + noise)
        plot(Function(Vh[STATE], u_o), title = "Observation")
    
    def cost(self, x):
        """
        Given the list x = [u,a,p] which describes the state, control, and
        adjoint variable compute the cost functional as the sum of 
        the misfit functional and the regularization functional.
        
        Return the list [cost functional, regularization functional, misfit functional]
        
        Note: p is not needed to compute the cost functional
        """        
        assert x[STATE] != None
                
        diff = x[STATE] - self.u_o
        Wuudiff = self.Wuu*diff
        misfit = .5 * diff.inner(Wuudiff)
        
        Rx = Vector()
        self.Prior.init_vector(Rx,0)
        self.Prior.R.mult(x[CONTROL], Rx)
        reg = .5 * x[CONTROL].inner(Rx)
        
        c = misfit + reg
        
        return c, reg, misfit
    
    def solveFwd(self, out, x, tol=1e-9):
        """
        Solve the forward problem.
        """
        A, b = self.assembleA(x, assemble_rhs = True)
        A.init_vector(out, 1)
        solver = PETScKrylovSolver("cg", "petsc_amg")
        solver.parameters["relative_tolerance"] = tol
        solver.set_operator(A)
        nit = solver.solve(out,b)
        
#        print "FWD", (self.A*out - b).norm("l2")/b.norm("l2"), nit

    
    def solveAdj(self, out, x, tol=1e-9):
        """
        Solve the adjoint problem.
        """
        At, badj = self.assembleA(x, assemble_adjoint = True,assemble_rhs = True)
        At.init_vector(out, 1)
        
        solver = PETScKrylovSolver("cg", "petsc_amg")
        solver.parameters["relative_tolerance"] = tol
        solver.set_operator(At)
        nit = solver.solve(out,badj)
        
#        print "ADJ", (self.At*out - badj).norm("l2")/badj.norm("l2"), nit
    
    def evalGradientControl(self,x, mg):
        """
        Evaluate the gradient for the variation control equation at the point x=[u,a,p].
        Parameters:
        - x = [u,a,p] the point at which to evaluate the gradient.
        - mg the variational gradient (g, atest) being atest a test function in the control space
          (Output parameter)
        
        Returns the norm of the gradient in the correct inner product g_norm = sqrt(g,g)
        """ 
        C = self.assembleC(x)

        self.Prior.init_vector(mg,0)
        C.transpmult(x[ADJOINT], mg)
        Rx = Vector()
        self.Prior.init_vector(Rx,0)
        self.Prior.R.mult(x[CONTROL], Rx)   
        mg.axpy(1., Rx)
        
        g = Vector()
        self.Prior.init_vector(g,1)
        
        solver = PETScKrylovSolver("cg", "jacobi")
        solver.set_operator(self.M)
        solver.solve(g, mg)
        g_norm = sqrt( g.inner(mg) )
        
        return g_norm
        
    
    def setPointForHessianEvaluations(self, x):  
        """
        Specify the point x = [u,a,p] at which the Hessian operator (or the Gauss-Newton approximation)
        need to be evaluated.
        """      
        self.A  = self.assembleA(x)
        self.At = self.assembleA(x, assemble_adjoint=True )
        self.C  = self.assembleC(x)
        self.Wau = self.assembleWau(x)
        self.Raa = self.assembleRaa(x)

        
    def solveFwdIncremental(self, sol, rhs, tol):
        """
        Solve the incremental forward problem for a given rhs
        """
        solver = PETScKrylovSolver("cg", "petsc_amg")
        solver.set_operator(self.A)
        solver.parameters["relative_tolerance"] = tol
        self.A.init_vector(sol,1)
        nit = solver.solve(sol,rhs)
#        print "FwdInc", (self.A*sol-rhs).norm("l2")/rhs.norm("l2"), nit
        
    def solveAdjIncremental(self, sol, rhs, tol):
        """
        Solve the incremental adjoint problem for a given rhs
        """
        solver = PETScKrylovSolver("cg", "petsc_amg")
        solver.set_operator(self.At)
        solver.parameters["relative_tolerance"] = tol
        self.At.init_vector(sol,1)
        nit = solver.solve(sol, rhs)
#        print "AdjInc", (self.At*sol-rhs).norm("l2")/rhs.norm("l2"), nit
    
    def applyC(self, da, out):
        self.C.mult(da,out)
    
    def applyCt(self, dp, out):
        self.C.transpmult(dp,out)
    
    def applyWuu(self, du, out):
        self.Wuu.mult(du, out)
    
    def applyWua(self, da, out):
        self.Wau.transpmult(da,out)

    
    def applyWau(self, du, out):
        self.Wau.mult(du, out)
    
    def applyR(self, da, out):
        self.Prior.R.mult(da, out)
        
    def Rsolver(self):        
        return self.Prior.Rsolver
    
    def applyRaa(self, da, out):
        self.Raa.mult(da, out)
            
if __name__ == "__main__":
    set_log_active(False)
    nx = 64
    ny = 64
    mesh = UnitSquareMesh(nx, ny)
    Vh2 = FunctionSpace(mesh, 'Lagrange', 2)
    Vh1 = FunctionSpace(mesh, 'Lagrange', 1)
    Vh = [Vh2, Vh1, Vh2]
    
    Prior = LaplacianPrior(Vh[CONTROL], gamma=1e-8, delta=1e-9)
    model = Poisson(mesh, Vh, Prior)
        
    a0 = interpolate(Expression("sin(x[0])"), Vh[CONTROL])
    modelVerify(model, a0.vector(), 1e-4, 1e-4)

    a0 = interpolate(Expression("0.0"),Vh[CONTROL])
    solver = ReducedSpaceNewtonCG(model)
    solver.parameters["abs_tolerance"] = 1e-9
    solver.parameters["inner_rel_tolerance"] = 1e-15
    solver.parameters["c_armijo"] = 1e-4
    solver.parameters["GN_iter"] = 5
    
    x = solver.solve(a0.vector())
    
    if solver.converged:
        print "\nConverged in ", solver.it, " iterations."
    else:
        print "\nNot Converged"

    print "Termination reason: ", solver.termination_reasons[solver.reason]
    print "Final gradient norm: ", solver.final_grad_norm
    print "Final cost: ", solver.final_cost
    
    xx = [Function(Vh[i], x[i]) for i in range(len(Vh))]
    plot(xx[STATE], title = "State")
    plot(exp(xx[CONTROL]), title = "exp(Control)")
    plot(xx[ADJOINT], title = "Adjoint")
    #interactive()
    
    model.setPointForHessianEvaluations(x)
    Hmisfit = ReducedHessian(model, solver.parameters["inner_rel_tolerance"], gauss_newton_approx=False, misfit_only=True)
    p = 50
    k = min( 250, Vh[CONTROL].dim()-p)
    Omega = np.random.randn(x[CONTROL].array().shape[0], k+p)
    d, U = singlePassG(Hmisfit, Prior.R, Prior.Rsolver, Omega, k)
    plt.plot(range(0,k), d, 'b*')
    plt.yscale('log')
    
    interactive()
    plt.show()
    
