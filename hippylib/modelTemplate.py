class ModelTemplate:
    """
    This class is a template for all the methods that a model object should
    provide.
    In the following we will denote with
    - u the state variable
    - a the control variable
    - p the adjoint variable
    
    For a concrete example see application/poisson/model.py.
    
    """
    
    def __init__(self, mesh, Vh, prior, other):
        """
        Construct a model by proving a mesh, the finite element spaces
        for the STATE/ADJOINT variable and the control variable, and a
        model for the prior information/regularization
        Pass any other parameter as needed.
        """
                
    def generate_vector(self, component = "ALL"):
        """
        By default, return the list [u,a,p] where:
        - u is any object that describes the state variable
        - a is a Vector object that describes the control variable.
          (Need to support linear algebra operations)
        - p is any object that describes the adjoint variable
        
        If component = STATE return only u
           component = CONTROL return only a
           component = ADJOINT return only p
        """ 
        return [None, None, None] #[u,a,p]
    
    def init_control(self, a):
        """
        Reshape a so that it is compatible with the control variable
        
        If R is a dolfin.Matrix object this can be simply implemented as
        model.R.init_vector(a,0)
        """
        return
            
    def cost(self, x):
        """
        Given the list x = [u,a,p] which describes the state, control, and
        adjoint variable compute the cost functional as the sum of 
        the misfit functional and the regularization functional.
        
        Return the list [cost functional, regularization functional, misfit functional]
        
        Note: p is not needed to compute the cost functional
        """
        return [None, None, None] #[cost, reg, misfit]
    
    def solveFwd(self, out, x, tol=1e-9):
        """
        Solve the (possibly non-linear) forward problem.
        Parameters:
        - out: is the solution of the forward problem (i.e. the state) (Output parameters)
        - x = [u,a,p] provides
              1) the control variable a for the solution of the forward problem
              2) the initial guess u if the forward problem is non-linear
          Note: p is not accessed
        - tol is the relative tolerance for the solution of the forward problem.
              [Default 1e-9].
        """
        return

    
    def solveAdj(self, out, x, tol=1e-9):
        """
        Solve the linear adjoint problem.
        Parameters:
        - out: is the solution of the adjoint problem (i.e. the adjoint p) (Output parameter)
        - x = [u,a,p] provides
              1) the control variable a for assembling the adjoint operator
              2) the state variable u for assembling the adjoint right hand side
          Note: p is not accessed
        - tol is the relative tolerance for the solution of the adjoint problem.
              [Default 1e-9].
        """
        return
    
    def evalGradientControl(self,x, mg):
        """
        Evaluate the gradient for the variation control equation at the point x=[u,a,p].
        Parameters:
        - x = [u,a,p] the point at which to evaluate the gradient.
        - mg the variational gradient (g, atest) being atest a test function in the control space
          (Output parameter)
        
        Returns the norm of the gradient in the correct inner product g_norm = sqrt(g,g)
        """ 
        return None #gradient norm
        
    
    def setPointForHessianEvaluations(self, x):
        """
        Specify the point x = [u,a,p] at which the Hessian operator (or the Gauss-Newton approximation)
        need to be evaluated.
        Parameters:
        - x = [u,a,p]: the point at which the Hessian or its Gauss-Newton approximation need to be evaluated.
        
        Note this routine should either:
        - simply store a copy of x and evaluate action of blocks of the Hessian on the fly
        - partially precompute the block of the hessian (if feasible)
        """
        return

        
    def solveFwdIncremental(self, sol, rhs, tol):
        """
        Solve the linearized (incremental) forward problem for a given rhs
        Parameters:
        - sol the solution of the linearized forward problem (Output)
        - rhs the right hand side of the linear system
        - tol the relative tolerance for the linear system
        """
        return
        
    def solveAdjIncremental(self, sol, rhs, tol):
        """
        Solve the incremental adjoint problem for a given rhs
        Parameters:
        - sol the solution of the incremental adjoint problem (Output)
        - rhs the right hand side of the linear system
        - tol the relative tolerance for the linear system
        """
        return
    
    def applyC(self, da, out):
        """
        Apply the C block of the Hessian to a (incremental) control variable.
        out = C da
        Parameters:
        - da the (incremental) control variable
        - out the action of the C block on da
        
        Note: this routine assumes that out has the correct shape.
        """
        return
    
    def applyCt(self, dp, out):
        """
        Apply the transpose of the C block of the Hessian to a (incremental) adjoint variable.
        out = C^t dp
        Parameters:
        - dp the (incremental) adjoint variable
        - out the action of the C^T block on dp
        
        Note: this routine assumes that out has the correct shape.
        """
        return

    
    def applyWuu(self, du, out):
        """
        Apply the Wuu block of the Hessian to a (incremental) state variable.
        out = Wuu du
        Parameters:
        - du the (incremental) state variable
        - out the action of the Wuu block on du
        
        Note: this routine assumes that out has the correct shape.
        """
        return
    
    def applyWua(self, da, out):
        """
        Apply the Wua block of the Hessian to a (incremental) control variable.
        out = Wua da
        Parameters:
        - da the (incremental) control variable
        - out the action of the Wua block on du
        
        Note: this routine assumes that out has the correct shape.
        """
        return

    
    def applyWau(self, du, out):
        """
        Apply the Wau block of the Hessian to a (incremental) state variable.
        out = Wau du
        Parameters:
        - du the (incremental) state variable
        - out the action of the Wau block on du
        
        Note: this routine assumes that out has the correct shape.
        """
        return
    
    def applyR(self, da, out):
        """
        Apply the regularization R to a (incremental) control variable.
        out = R da
        Parameters:
        - da the (incremental) control variable
        - out the action of R on da
        
        Note: this routine assumes that out has the correct shape.
        """
        return
    
    def RPreconditioner(self):
        """
        Return an object Rprec that is a suitable preconditioner for the regularization
        operator R.
        
        The preconditioner object should implement the method Rprec.solve(z,r) such that
        R*z \\approx r.
        """

    
    def applyRaa(self, da, out):
        """
        Apply the Raa block of the Hessian to a (incremental) control variable.
        out = Raa da
        Parameters:
        - da the (incremental) control variable
        - out the action of R on da
        
        Note: this routine assumes that out has the correct shape.
        """
        return
