import dolfin as dl
import numpy as np
from linalg import MatPtAP, MatMatMult, MatAtB, get_diagonal, estimate_diagonal_inv_coloring, getColoring, to_dense, amg_method
from traceEstimator import TraceEstimator
from assemblePointwiseObservation import assemblePointwiseObservation
import math
from expression import code_Mollifier
from cgsampler import CGSampler


class _RinvM:
    """
    Operator that models the action of R^{-1}M.
    It is used in the randomized trace estimator.
    """
    def __init__(self, Rsolver, M):
        self.Rsolver = Rsolver
        self.M = M
        
    def init_vector(self,x,dim):
        self.M.init_vector(x,dim)
            
    def mult(self,x,y):
        self.Rsolver.solve(y, self.M*x)

class _Prior:
    """
    Abstract class to describe the prior model.
    Concrete instances of a _Prior class should expose
    the following attributes and methods.
    
    attributes:
    - R:       an operator to apply the regularization/precision
               operator
    - Rsolver: an operator to apply the inverse of the
               regularization/precision operator
    - M:       the mass matrix in the control space
    - mean:    the prior mean
    
    methods:
    - init_vector(self,x,dim):
      Inizialize a vector x to be compatible with the range/domain of R
      If dim == "noise" inizialize x to be compatible with the size of
      white noise used for sampling.
      
    - sample(self, noise, s, add_mean=True):
      Given a noise ~ N(0, I) compute a sample s from the prior.
      If add_mean=True add the prior mean value to s.
    """ 
               
    def trace(self, method="Exact", tol=1e-1, min_iter=20, max_iter=100):
        """
        Compute/Estimate the trace of the prior covariance operator.
        
        - If method=="Exact" we compute the trace exactly by summing the diagonal entries of R^{-1}M.
          This requires to solve n linear system in R (not scalable, but ok for illustration purposes).
          
        - If method=="Estimator" use the trace estimator algorithms implemeted in the class TraceEstimator.
          tol is a relative bound on the estimator standard deviation. In particular, we used enough samples in the
          Estimator such the std deviation of the estimator is less then tol*tr(Prior).
          min_iter and max_iter are the lower and upper bound on the number of samples to be used for the
          estimation of the trace. 
        """
        op = _RinvM(self.Rsolver, self.M)
        if method == "Exact":
            marginal_variance = dl.Vector()
            self.init_vector(marginal_variance,0)
            get_diagonal(op, marginal_variance, solve_mode=False)
            return marginal_variance.sum()
        elif method == "Estimator":
            tr_estimator = TraceEstimator(op, False, tol)
            tr_exp, tr_var = tr_estimator(min_iter, max_iter)
            return tr_exp
        else:
            raise NameError("Unknown method")
        
    def pointwise_variance(self, method, path_len = 8):
        """
        Compute/Estimate the prior pointwise variance.
        
        - If method=="Exact" we compute the diagonal entries of R^{-1} entry by entry. 
          This requires to solve n linear system in R (not scalable, but ok for illustration purposes).
          
        - If method=="ProbingEstimator" we use a probing algorithm to approximate the diagonal entries of R^{-1}.
          The number of linear system in R to solve depends on the sparsity pattern of R.
          Even if the computational cost is much lower than the exact method, the complexity still increase linearly
          with the size of the problem. path_len represent the power of R used to compute the adjacency matrix for the
          graph coloring algorithm. See function estimate_diagonal_inv_coloring for details.
        """
        pw_var = dl.Vector()
        self.init_vector(pw_var,0)
        if method == "Exact":
            get_diagonal(self.Rsolver, pw_var, solve_mode=True)
        elif method == "ProbingEstimator":
            if type(self.R) == dl.Matrix:
                coloring = getColoring(self.R, path_len)
            else:
                Rsparsity = MatPtAP(self.M, self.A)
                coloring = getColoring(Rsparsity, path_len)
                
            estimate_diagonal_inv_coloring(self.Rsolver, coloring, pw_var)
        else:
            raise NameError("Unknown method")
        
        return pw_var
    
    def cost(self,a):
        d = self.mean.copy()
        d.axpy(-1., a)
        Rd = dl.Vector()
        self.init_vector(Rd,0)
        self.R.mult(d,Rd)
        return .5*Rd.inner(d)
    
    def grad(self,a, out):
        d = a.copy()
        d.axpy(-1., self.mean)
        self.R.mult(d,out)

class LaplacianPrior(_Prior):
    """
    This class implement a Prior model with covariance matrix
    C = (\delta I + \gamma \Delta) ^ {-1}.
    
    The magnitude of \gamma governs the variance of the samples, while
    the ratio \frac{\gamma}{\delta} governs the correlation lenght.
    
    Note that C is a trace class operator only in 1D while it is not
    a valid prior in 2D and 3D.
    """
    
    def __init__(self, Vh, gamma, delta, mean=None, rel_tol=1e-12, max_iter=100):
        """
        Construct the Prior model.
        Input:
        - Vh:              the finite element space for the parameter
        - gamma and delta: the coefficient in the PDE
        - Theta:           the s.p.d. tensor for anisotropic diffusion of the pde
        - mean:            the prior mean
        """        
        assert delta != 0., "Intrinsic Gaussian Prior are not supported"
        self.Vh = Vh
        
        trial = dl.TrialFunction(Vh)
        test  = dl.TestFunction(Vh)
        
        varfL = dl.inner(dl.nabla_grad(trial), dl.nabla_grad(test))*dl.dx
        varfM = dl.inner(trial,test)*dl.dx
        
        self.M = dl.assemble(varfM)
        self.R = dl.assemble(gamma*varfL + delta*varfM)
        
        self.Rsolver = dl.PETScKrylovSolver("cg", amg_method())
        self.Rsolver.set_operator(self.R)
        self.Rsolver.parameters["maximum_iterations"] = max_iter
        self.Rsolver.parameters["relative_tolerance"] = rel_tol
        self.Rsolver.parameters["error_on_nonconvergence"] = True
        self.Rsolver.parameters["nonzero_initial_guess"] = False
        
        self.Msolver = dl.PETScKrylovSolver("cg", "jacobi")
        self.Msolver.set_operator(self.M)
        self.Msolver.parameters["maximum_iterations"] = max_iter
        self.Msolver.parameters["relative_tolerance"] = rel_tol
        self.Msolver.parameters["error_on_nonconvergence"] = True
        self.Msolver.parameters["nonzero_initial_guess"] = False
        
        Q1h = dl.FunctionSpace(Vh.mesh(), 'Quadrature', 2*Vh._FunctionSpace___degree)
        ndim = Vh.mesh().geometry().dim()
        Qh = dl.MixedFunctionSpace([Q1h for i in range(ndim+1)])
        ph = dl.TrialFunction(Qh)
        qh = dl.TestFunction(Qh)
        
        pph = dl.split(ph)
        
        Mqh = dl.assemble(dl.inner(ph, qh)*dl.dx)
        ones = dl.Vector()
        Mqh.init_vector(ones,0)
        ones.set_local( np.ones(ones.array().shape, dtype =ones.array().dtype ) )
        dMqh = Mqh*ones
        dMqh.set_local( ones.array() / np.sqrt(dMqh.array() ) )
        Mqh.zero()
        Mqh.set_diagonal(dMqh)
        
        sqrtdelta = math.sqrt(delta)
        sqrtgamma = math.sqrt(gamma)
        varfGG = sqrtdelta*pph[0]*test*dl.dx
        for i in range(ndim):
            varfGG = varfGG + sqrtgamma*pph[i+1]*test.dx(i)*dl.dx
            
        GG = dl.assemble(varfGG)
        self.sqrtR = MatMatMult(GG, Mqh)
                        
        self.mean = mean
        
        if self.mean is None:
            self.mean = dl.Vector()
            self.init_vector(self.mean, 0)
        
    def init_vector(self,x,dim):
        """
        Inizialize a vector x to be compatible with the range/domain of R.
        If dim == "noise" inizialize x to be compatible with the size of
        white noise used for sampling.
        """
        if dim == "noise":
            self.sqrtR.init_vector(x,1)
        else:
            self.R.init_vector(x,dim)
                
    def sample(self, noise, s, add_mean=True):
        """
        Given a noise ~ N(0, I) compute a sample s from the prior.
        If add_mean=True add the prior mean value to s.
        """

        rhs = self.sqrtR*noise
        self.Rsolver.solve(s,rhs)
        
        if add_mean:
            s.axpy(1., self.mean)
        

class _BilaplacianR:
    """
    Operator that represent the action of the regularization/precision matrix
    for the Bilaplacian prior.
    """
    def __init__(self, A, Msolver):
        self.A = A
        self.Msolver = Msolver

        self.help1, self.help2 = dl.Vector(), dl.Vector()
        self.A.init_vector(self.help1, 0)
        self.A.init_vector(self.help2, 1)
        
    def init_vector(self,x, dim):
        self.A.init_vector(x,1)
        
    def inner(self,x,y):
        Rx = dl.Vector()
        self.init_vector(Rx,0)
        self.mult(x, Rx)
        return Rx.inner(y)
        
    def mult(self,x,y):
        self.A.mult(x, self.help1)
        self.Msolver.solve(self.help2, self.help1)
        self.A.mult(self.help2, y)
        
class _BilaplacianRsolver():
    """
    Operator that represent the action of the inverse the regularization/precision matrix
    for the Bilaplacian prior.
    """
    def __init__(self, Asolver, M):
        self.Asolver = Asolver
        self.M = M
        
        self.help1, self.help2 = dl.Vector(), dl.Vector()
        self.init_vector(self.help1, 0)
        self.init_vector(self.help2, 0)
        
    def init_vector(self,x, dim):
        self.M.init_vector(x,1)
        
    def solve(self,x,b):
        nit = self.Asolver.solve(self.help1, b)
        self.M.mult(self.help1, self.help2)
        nit += self.Asolver.solve(x, self.help2)
        return nit
        
        
class BiLaplacianPrior(_Prior):
    """
    This class implement a Prior model with covariance matrix
    C = (\delta I + \gamma \div Theta \grad) ^ {-2}.
    
    The magnitude of \delta\gamma governs the variance of the samples, while
    the ratio \frac{\gamma}{\delta} governs the correlation lenght.
    
    Here Theta is a s.p.d tensor that models anisotropy in the covariance kernel.
    """
    
    def __init__(self, Vh, gamma, delta, Theta = None, mean=None, rel_tol=1e-12, max_iter=1000):
        """
        Construct the Prior model.
        Input:
        - Vh:              the finite element space for the parameter
        - gamma and delta: the coefficient in the PDE
        - Theta:           the s.p.d. tensor for anisotropic diffusion of the pde
        - mean:            the prior mean
        """
        assert delta != 0., "Intrinsic Gaussian Prior are not supported"
        self.Vh = Vh
        
        trial = dl.TrialFunction(Vh)
        test  = dl.TestFunction(Vh)
        
        if Theta == None:
            varfL = dl.inner(dl.nabla_grad(trial), dl.nabla_grad(test))*dl.dx
        else:
            varfL = dl.inner( Theta*dl.grad(trial), dl.grad(test))*dl.dx
        
        varfM = dl.inner(trial,test)*dl.dx
        
        self.M = dl.assemble(varfM)
        self.Msolver = dl.PETScKrylovSolver("cg", "jacobi")
        self.Msolver.set_operator(self.M)
        self.Msolver.parameters["maximum_iterations"] = max_iter
        self.Msolver.parameters["relative_tolerance"] = rel_tol
        self.Msolver.parameters["error_on_nonconvergence"] = True
        self.Msolver.parameters["nonzero_initial_guess"] = False
        
        self.A = dl.assemble(gamma*varfL + delta*varfM)        
        self.Asolver = dl.PETScKrylovSolver("cg", amg_method())
        self.Asolver.set_operator(self.A)
        self.Asolver.parameters["maximum_iterations"] = max_iter
        self.Asolver.parameters["relative_tolerance"] = rel_tol
        self.Asolver.parameters["error_on_nonconvergence"] = True
        self.Asolver.parameters["nonzero_initial_guess"] = False
        
        Qh = dl.FunctionSpace(Vh.mesh(), 'Quadrature', 2*Vh._FunctionSpace___degree)
        ph = dl.TrialFunction(Qh)
        qh = dl.TestFunction(Qh)
        Mqh = dl.assemble(ph*qh*dl.dx)
        ones = dl.interpolate(dl.Constant(1.), Qh).vector()
        dMqh = Mqh*ones
        Mqh.zero()
        dMqh.set_local( ones.array() / np.sqrt(dMqh.array() ) )
        Mqh.set_diagonal(dMqh)
        MixedM = dl.assemble(ph*test*dl.dx)
        self.sqrtM = MatMatMult(MixedM, Mqh)
                     
        self.R = _BilaplacianR(self.A, self.Msolver)      
        self.Rsolver = _BilaplacianRsolver(self.Asolver, self.M)
         
        self.mean = mean
        
        if self.mean is None:
            self.mean = dl.Vector()
            self.init_vector(self.mean, 0)
     
    def init_vector(self,x,dim):
        """
        Inizialize a vector x to be compatible with the range/domain of R.
        If dim == "noise" inizialize x to be compatible with the size of
        white noise used for sampling.
        """
        if dim == "noise":
            self.sqrtM.init_vector(x, 1)
        else:
            self.A.init_vector(x,dim)   
        
    def sample(self, noise, s, add_mean=True):
        """
        Given a noise ~ N(0, I) compute a sample s from the prior.
        If add_mean=True add the prior mean value to s.
        """
        rhs = self.sqrtM*noise
        self.Asolver.solve(s, rhs)
        
        if add_mean:
            s.axpy(1., self.mean)
        
class ConstrainedBiLaplacianPrior(_Prior):
    """
    This class implement a Prior model with covariance matrix
    C = (\delta I + \gamma \div Theta \grad + \sum_i dirac(x - x_i) ) ^ {-2}.
    
    Here x_i (i=1,...,n) are points were we assume to know exactly the value
    of the parameter (i.e. m(x_i) = m_true( x_i) for i=1,...,n).
    
    The magnitude of \delta\gamma governs the variance of the samples, while
    the ratio \frac{\gamma}{\delta} governs the correlation lenght.
    
    Here Theta is a s.p.d tensor that models anisotropy in the covariance kernel.
    
    The prior mean is computed by solving the least square problem
    min || m ||^2_R, s.t. m(x_i) = m_true(x_i) for (i=1,...,n).
    """
    
    def __init__(self, Vh, gamma, delta, locations, m_true, Theta = None, pen = 1e4, rel_tol=1e-12, max_iter=1000):
        """
        Construct the Prior model.
        Input:
        - Vh:              the finite element space for the parameter
        - gamma and delta: the coefficient in the PDE
        - locations:       the points x_i at which we assume to know the
                           true value of the parameter
        - m_true:          the true model
        - Theta:           the s.p.d. tensor for anisotropic diffusion of the pde
        - pen:             a discretization parameter for the dirac delta
        """
        assert delta != 0. or pen != 0, "Intrinsic Gaussian Prior are not supported"
        self.Vh = Vh
        
        trial = dl.TrialFunction(Vh)
        test  = dl.TestFunction(Vh)
        
        if Theta == None:
            varfL = dl.inner(dl.nabla_grad(trial), dl.nabla_grad(test))*dl.dx
        else:
            varfL = dl.inner(Theta*dl.grad(trial), dl.grad(test))*dl.dx
        varfM = dl.inner(trial,test)*dl.dx
        
        self.M = dl.assemble(varfM)
        self.Msolver = dl.PETScKrylovSolver("cg", "jacobi")
        self.Msolver.set_operator(self.M)
        self.Msolver.parameters["maximum_iterations"] = max_iter
        self.Msolver.parameters["relative_tolerance"] = rel_tol
        self.Msolver.parameters["error_on_nonconvergence"] = True
        self.Msolver.parameters["nonzero_initial_guess"] = False
        
        Po = assemblePointwiseObservation(Vh, locations)
                
        self.A = dl.assemble(gamma*varfL+delta*varfM)
        PotPo = MatAtB(Po,Po)
        self.A.axpy(pen, PotPo, False);      
        self.Asolver = dl.PETScKrylovSolver("cg", amg_method())
        self.Asolver.set_operator(self.A)
        self.Asolver.parameters["maximum_iterations"] = max_iter
        self.Asolver.parameters["relative_tolerance"] = rel_tol
        self.Asolver.parameters["error_on_nonconvergence"] = True
        self.Asolver.parameters["nonzero_initial_guess"] = False
        
        Qh = dl.FunctionSpace(Vh.mesh(), 'Quadrature', 2*Vh._FunctionSpace___degree)
        ph = dl.TrialFunction(Qh)
        qh = dl.TestFunction(Qh)
        Mqh = dl.assemble(ph*qh*dl.dx)
        ones = dl.interpolate(dl.Constant(1.), Qh).vector()
        dMqh = Mqh*ones
        Mqh.zero()
        dMqh.set_local( ones.array() / np.sqrt(dMqh.array() ) )
        Mqh.set_diagonal(dMqh)
        MixedM = dl.assemble(ph*test*dl.dx)
        self.sqrtM = MatMatMult(MixedM, Mqh)
             
        self.R = _BilaplacianR(self.A, self.Msolver)      
        self.Rsolver = _BilaplacianRsolver(self.Asolver, self.M)
        
        rhs = dl.Vector()
        self.mean = dl.Vector()
        self.init_vector(rhs, 0)
        self.init_vector(self.mean, 0)
        
        PotPo.mult(m_true, rhs)
        rhs *= pen  
        self.Asolver.solve(self.mean, rhs)
        
     
    def init_vector(self,x,dim):
        """
        Inizialize a vector x to be compatible with the range/domain of R.
        If dim == "noise" inizialize x to be compatible with the size of
        white noise used for sampling.
        """
        if dim == "noise":
            self.sqrtM.init_vector(x, 1)
        else:
            self.A.init_vector(x,dim)   
        
    def sample(self, noise, s, add_mean=True):
        """
        Given a noise ~ N(0, I) compute a sample s from the prior.
        If add_mean=True add the prior mean value to s.
        """
        rhs = self.sqrtM*noise
        self.Asolver.solve(s, rhs)
        
        if add_mean:
            s.axpy(1., self.mean)

class MollifiedBiLaplacianPrior(_Prior):
    """
    This class implement a Prior model with covariance matrix
    C = ( (\delta + pen \sum_i m(x - x_i) ) I + \gamma \div Theta \grad) ^ {-2},
    
    where
    - Theta is a s.p.d tensor that models anisotropy in the covariance kernel
    - x_i (i=1,...,n) are points were we assume to know exactly the value
    of the parameter (i.e. m(x_i) = m_true( x_i) for i=1,...,n).    
    - m is the mollifier function: m(x - x_i) = exp( - [\frac{\gamma}{\delta}|| x - x_i ||_Theta^{-1}]^order ).
    - pen is a penalization paramter
    
    The magnitude of \delta\gamma governs the variance of the samples, while
    the ratio \frac{\gamma}{\delta} governs the correlation lenght.
    
    The prior mean is computed by solving 
    ( (\delta + \sum_i m(x - x_i) ) I + \gamma \div Theta \grad ) m = \sum_i m(x - x_i) m_true.
    """
    
    def __init__(self, Vh, gamma, delta, locations, m_true, Theta = None, pen = 1e1, order=2, rel_tol=1e-12, max_iter=1000):
        """
        Construct the Prior model.
        Input:
        - Vh:              the finite element space for the parameter
        - gamma and delta: the coefficient in the PDE
        - locations:       the points x_i at which we assume to know the
                           true value of the parameter
        - m_true:          the true model
        - Theta:           the s.p.d. tensor for anisotropic diffusion of the pde
        - pen:             a penalization parameter for the mollifier
        """
        assert delta != 0. or pen != 0, "Intrinsic Gaussian Prior are not supported"
        self.Vh = Vh
        
        trial = dl.TrialFunction(Vh)
        test  = dl.TestFunction(Vh)
        
        if Theta == None:
            varfL = dl.inner(dl.nabla_grad(trial), dl.nabla_grad(test))*dl.dx
        else:
            varfL = dl.inner(Theta*dl.grad(trial), dl.grad(test))*dl.dx
        varfM = dl.inner(trial,test)*dl.dx
        
        self.M = dl.assemble(varfM)
        self.Msolver = dl.PETScKrylovSolver("cg", "jacobi")
        self.Msolver.set_operator(self.M)
        self.Msolver.parameters["maximum_iterations"] = max_iter
        self.Msolver.parameters["relative_tolerance"] = rel_tol
        self.Msolver.parameters["error_on_nonconvergence"] = True
        self.Msolver.parameters["nonzero_initial_guess"] = False
        
        #mfun = Mollifier(gamma/delta, dl.inv(Theta), order, locations)
        mfun = dl.Expression(code_Mollifier)
        mfun.l = gamma/delta
        mfun.o = order
        mfun.theta0 = 1./Theta.theta0
        mfun.theta1 = 1./Theta.theta1
        mfun.alpha = Theta.alpha
        for ii in range(locations.shape[0]):
            mfun.addLocation(locations[ii,0], locations[ii,1])
            
        varfmo = mfun*dl.inner(trial,test)*dl.dx
        MO = dl.assemble(pen*varfmo)
                
        self.A = dl.assemble(gamma*varfL+delta*varfM + pen*varfmo)
     
        self.Asolver = dl.PETScKrylovSolver("cg", amg_method())
        self.Asolver.set_operator(self.A)
        self.Asolver.parameters["maximum_iterations"] = max_iter
        self.Asolver.parameters["relative_tolerance"] = rel_tol
        self.Asolver.parameters["error_on_nonconvergence"] = True
        self.Asolver.parameters["nonzero_initial_guess"] = False
        
        Qh = dl.FunctionSpace(Vh.mesh(), 'Quadrature', 2*Vh._FunctionSpace___degree)
        ph = dl.TrialFunction(Qh)
        qh = dl.TestFunction(Qh)
        Mqh = dl.assemble(ph*qh*dl.dx)
        ones = dl.interpolate(dl.Constant(1.), Qh).vector()
        dMqh = Mqh*ones
        Mqh.zero()
        dMqh.set_local( ones.array() / np.sqrt(dMqh.array() ) )
        Mqh.set_diagonal(dMqh)
        MixedM = dl.assemble(ph*test*dl.dx)
        self.sqrtM = MatMatMult(MixedM, Mqh)
             
        self.R = _BilaplacianR(self.A, self.Msolver)      
        self.Rsolver = _BilaplacianRsolver(self.Asolver, self.M)
        
        rhs = dl.Vector()
        self.mean = dl.Vector()
        self.init_vector(rhs, 0)
        self.init_vector(self.mean, 0)
        
        MO.mult(m_true, rhs)
        self.Asolver.solve(self.mean, rhs)
        
     
    def init_vector(self,x,dim):
        """
        Inizialize a vector x to be compatible with the range/domain of R.
        If dim == "noise" inizialize x to be compatible with the size of
        white noise used for sampling.
        """
        if dim == "noise":
            self.sqrtM.init_vector(x, 1)
        else:
            self.A.init_vector(x,dim)   
        
    def sample(self, noise, s, add_mean=True):
        """
        Given a noise ~ N(0, I) compute a sample s from the prior.
        If add_mean=True add the prior mean value to s.
        """
        rhs = self.sqrtM*noise
        self.Asolver.solve(s, rhs)
        
        if add_mean:
            s.axpy(1., self.mean)


class LaplaceBeltramiPrior(_Prior):
    def __init__(self, Vh, gamma, delta, ds, mean=None, rel_tol=1e-12, max_iter=100):
        
        assert delta != 0., "Intrinsic Gaussian Prior are not supported"
        self.Vh = Vh
        trial = dl.TrialFunction(Vh)
        test  = dl.TestFunction(Vh)
        
        n = dl.FacetNormal(Vh.mesh()) 
        def dir_grad(uh,n):
            return dl.grad(uh) - dl.dot(dl.outer(n,n), dl.grad(uh))
        
        varfL = gamma*dl.inner(dir_grad(trial,n), dir_grad(test,n))*ds
        varfM = delta*dl.inner(trial,test)*ds
        
        self.M = dl.assemble(varfM)
        self.R = dl.assemble(varfL + varfM)
            
        self.M_1 = dl.assemble(varfM, keep_diagonal=True)
        self.R_1 = dl.assemble(varfL + varfM, keep_diagonal=True)
        self.M_1.ident_zeros()
        self.R_1.ident_zeros()
            
        self.Rsolver = dl.PETScKrylovSolver("cg", amg_method())
        self.Rsolver.set_operator(self.R_1)
        self.Rsolver.parameters["maximum_iterations"] = max_iter
        self.Rsolver.parameters["relative_tolerance"] = rel_tol
        self.Rsolver.parameters["error_on_nonconvergence"] = True
        self.Rsolver.parameters["nonzero_initial_guess"] = False
        
        self.Msolver = dl.PETScKrylovSolver("cg", "jacobi")
        self.Msolver.set_operator(self.M_1)
        self.Msolver.parameters["maximum_iterations"] = max_iter
        self.Msolver.parameters["relative_tolerance"] = rel_tol
        self.Msolver.parameters["error_on_nonconvergence"] = True
        self.Msolver.parameters["nonzero_initial_guess"] = False
        
        self.sampler = CGSampler()
        self.sampler.set_operator(self.R_1)
        
        self.mean = mean
        
        if self.mean is None:
            self.mean = dl.Vector()
            self.init_vector(self.mean, 0)
    
    def init_vector(self,x,dim):
        """
        Inizialize a vector x to be compatible with the range/domain of R.
        If dim == "noise" inizialize x to be compatible with the size of
        white noise used for sampling.
        """
        if dim == "noise":
            self.R.init_vector(x, 1)
        else:
            self.R.init_vector(x,dim)
            
    def sample(self, noise, s, add_mean=True):
        """
        Given a noise ~ N(0, I) compute a sample s from the prior.
        If add_mean=True add the prior mean value to s.
        """
        self.sampler.sample(noise.array(), s)
        
        if add_mean:
            s.axpy(1., self.mean)         