from variables import STATE, PARAMETER, ADJOINT
from linalg import MatPtAP
from dolfin import Vector, PETScKrylovSolver

class ReducedHessian:
    """
    This class implements matrix free application of the reduced hessian operator.
    The constructor takes the following parameters:
    - model:               the object which contains the description of the problem.
    - innerTol:            the relative tolerance for the solution of the incremental
                           forward and adjoint problems.
    - gauss_newton_approx: a boolean flag that describes whenever the true hessian or
                           the Gauss Newton approximation of the Hessian should be
                           applied.
    - misfit_only:         a boolean flag that describes whenever the full hessian or
                           only the misfit component of the hessian is used.
    
    Type help(modelTemplate) for more information on which methods model should implement.
    """
    def __init__(self, model, innerTol, gauss_newton_approx=False, misfit_only=False):
        """
        Construct the reduced Hessian Operator:
        """
        self.model = model
        self.tol = innerTol
        self.gauss_newton_approx = gauss_newton_approx
        self.misfit_only=misfit_only
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
        Apply the reduced Hessian (or the Gauss Newton approximation) to the vector x
        Return the result in y.
        """
        if self.gauss_newton_approx:
            self.GNHessian(x,y)
        else:
            self.TrueHessian(x,y)
        
        self.ncalls += 1
    
    def inner(self,x,y):
        """
        Perform the inner product between x and y in the norm induced by the reduced
        Hessian H.
        (x, y)_H = x' H y
        """
        Ay = self.model.generate_vector(PARAMETER)
        Ay.zero()
        self.mult(y,Ay)
        return x.inner(Ay)
            
    def GNHessian(self,x,y):
        """
        Apply the Gauss Newton approximation of the reduced Hessian to the vector x
        Return the result in y.        
        """
        self.model.applyC(x, self.rhs_fwd)
        self.model.solveFwdIncremental(self.uhat, self.rhs_fwd, self.tol)
        self.model.applyWuu(self.uhat, self.rhs_adj, True)
        self.model.solveAdjIncremental(self.phat, self.rhs_adj, self.tol)
        self.model.applyCt(self.phat, y)
        
        if not self.misfit_only:
            self.model.applyR(x,self.yhelp)
            y.axpy(1., self.yhelp)

        
    def TrueHessian(self, x, y):
        """
        Apply the the reduced Hessian to the vector x.
        Return the result in y.        
        """
        self.model.applyC(x, self.rhs_fwd)
        self.model.solveFwdIncremental(self.uhat, self.rhs_fwd, self.tol)
        self.model.applyWuu(self.uhat, self.rhs_adj, False)
        self.model.applyWua(x, self.rhs_adj2)
        self.rhs_adj.axpy(-1., self.rhs_adj2)
        self.model.solveAdjIncremental(self.phat, self.rhs_adj, self.tol)
        self.model.applyRaa(x, y)
        self.model.applyCt(self.phat, self.yhelp)
        y.axpy(1., self.yhelp)
        self.model.applyWau(self.uhat, self.yhelp)
        y.axpy(-1., self.yhelp)
        
        if not self.misfit_only:
            self.model.applyR(x,self.yhelp)
            y.axpy(1., self.yhelp)
        
class ReducedHessianActiveSet:
    """
    Experimental implementation of ActiveSet constraints
    """
    def __init__(self, model, innerTol, gauss_newton_approx=False):
        self.H = ReducedHessian(model, innerTol, gauss_newton_approx)
        self.I_set = model.getIdentityMatrix(PARAMETER)
        self.A_set = model.getIdentityMatrix(PARAMETER)
        self.A_set.zero()
        
        self.index_active_set = None
        self.x_inactive = Vector()
        self.y_inactive = Vector()
        self.I_diag = Vector()
        
        self.I_set.init_vector(self.x_inactive, 0)
        self.I_set.init_vector(self.y_inactive, 0)
        self.I_set.init_vector(self.I_diag, 0)
        self.active_set = False
        
    def init_vector(self, x, dim):
        self.H.init_vector(x, dim)
        
    def setActiveSet(self, index_active_set):
        self.index_active_set = index_active_set
        
        self.I_diag.set_local(abs( index_active_set ))
                       
        self.A_set.set_diagonal(self.I_diag)
        self.I_set.set_diagonal(1. - self.I_diag )
        
        self.active_set = True
        
    def applyBounds(self, x, lb, ub):
        # Only works for lb = 0, ub = None 
        if self.index_active_set is None:
            return
        
        if lb is not None:
            index = self.index_active_set == 1.
            x.array()[index] = lb[index] 
        if ub is not None:
            index = self.index_active_set == -1.
            x.array()[index] = ub[index] 
            
    def getPreconditioner(self):
        
        if self.active_set == False:
            return self.H.model.RPreconditioner()
        
        R_inactive = MatPtAP(self.H.model.Prior.R, self.I_set)
        R_inactive.ident_zeros()
#        R_inactive.axpy(1., self.A_set, False)
        
        Rprec = PETScKrylovSolver("richardson", "amg")
        Rprec.parameters["maximum_iterations"] = 1
        Rprec.parameters["error_on_nonconvergence"] = False
        Rprec.parameters["nonzero_initial_guess"] = False
        Rprec.set_operator(R_inactive)
        
        return Rprec
              
    def mult(self,x,y):
        if self.active_set:
            # y(active_set) = x(active_set)
            self.A_set.mult(x,y)
            # y(inactive_set)= H(inactive_set,inactive_set)*x(inactive_set)
            self.I_set.mult(x, self.x_inactive)
            self.H.mult(self.x_inactive, self.y_inactive)
            self.I_set.mult(self.y_inactive, self.x_inactive)
            y.axpy(1., self.x_inactive)
        else:                  
            self.H.mult(x, y)   
            
