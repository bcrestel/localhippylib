import dolfin as dl
from cgsampler import CGSampler
from matFreeChol import MatFreeChol
import numpy as np
from linalg import MatPtAP, MatMatMult, MatAtB, get_diagonal, estimate_diagonal_inv_coloring, getColoring, to_dense
from traceEstimator import TraceEstimator
from assemblePointwiseObservation import assemblePointwiseObservation
import math
from expression import code_AnisTensor2D, code_Mollifier


class _RinvM:
    def __init__(self, Rsolver, M):
        self.Rsolver = Rsolver
        self.M = M
        
    def init_vector(self,x,dim):
        self.M.init_vector(x,dim)
            
    def mult(self,x,y):
        self.Rsolver.solve(y, self.M*x)

class _Prior:            
    def trace(self, method="Exact", tol=1e-1, min_iter=20, max_iter=100):
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
    

class LaplacianPrior(_Prior):
    
    def __init__(self, Vh, gamma, delta, mean=None, rel_tol=1e-12, max_iter=100):
        
        assert delta != 0., "Intrinsic Gaussian Prior are not supported"
        self.Vh = Vh
        
        trial = dl.TrialFunction(Vh)
        test  = dl.TestFunction(Vh)
        
        varfL = gamma*dl.inner(dl.nabla_grad(trial), dl.nabla_grad(test))*dl.dx
        varfM = delta*dl.inner(trial,test)*dl.dx
        
        self.M = dl.assemble(varfM)
        self.R = dl.assemble(varfL + varfM)
        
        self.Rsolver = dl.PETScKrylovSolver("cg", "amg")
        self.Rsolver.set_operator(self.R)
        self.Rsolver.parameters["maximum_iterations"] = max_iter
        self.Rsolver.parameters["relative_tolerance"] = rel_tol
        self.Rsolver.parameters["error_on_nonconvergence"] = True
        self.Rsolver.parameters["nonzero_initial_guess"] = False
        
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
        if dim == "noise":
            self.sqrtR.init_vector(x,1)
        else:
            self.R.init_vector(x,dim)
                
    def sample(self, noise, s, add_mean=True):
        #sampler = CGSampler()
        #sampler.set_operator(self.R)
        #sampler.sample(noise.array(), s)
        rhs = self.sqrtR*noise
        self.Rsolver.solve(s,rhs)
        
        if add_mean:
            s.axpy(1., self.mean)
        

class _BilaplacianR:
    def __init__(self, A, Msolver):
        self.A = A
        self.Msolver = Msolver

        self.help1, self.help2 = dl.Vector(), dl.Vector()
        self.A.init_vector(self.help1, 0)
        self.A.init_vector(self.help2, 1)
        
    def init_vector(self,x, dim):
        self.A.init_vector(x,1)
        
    def mult(self,x,y):
        self.A.mult(x, self.help1)
        self.Msolver.solve(self.help2, self.help1)
        self.A.mult(self.help2, y)
        
class _BilaplacianRsolver():
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
    
    def __init__(self, Vh, gamma, delta, anis_diff = None, mean=None, rel_tol=1e-12, max_iter=1000):
        
        assert delta != 0., "Intrinsic Gaussian Prior are not supported"
        self.Vh = Vh
        
        trial = dl.TrialFunction(Vh)
        test  = dl.TestFunction(Vh)
        
        if anis_diff == None:
            varfL = dl.inner(dl.nabla_grad(trial), dl.nabla_grad(test))*dl.dx
        else:
            varfL = dl.inner( anis_diff*dl.grad(trial), dl.grad(test))*dl.dx
        
        varfM = dl.inner(trial,test)*dl.dx
        
        self.M = dl.assemble(varfM)
        self.Msolver = dl.PETScKrylovSolver("cg", "jacobi")
        self.Msolver.set_operator(self.M)
        self.Msolver.parameters["maximum_iterations"] = max_iter
        self.Msolver.parameters["relative_tolerance"] = rel_tol
        self.Msolver.parameters["error_on_nonconvergence"] = True
        self.Msolver.parameters["nonzero_initial_guess"] = False
        
        self.A = dl.assemble(gamma*varfL + delta*varfM)        
        self.Asolver = dl.PETScKrylovSolver("cg", "amg")
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
        
        
        #ones = dl.Vector()
        #self.M.init_vector(ones,0)
        #ones.set_local(np.ones(self.M.size(0)))
        #Mones = self.M*ones

        #self.sqrtM = dl.assemble(varfM)
        #s = dl.Vector()
        #self.M.init_vector(s,0)
        #s.set_local( np.sqrt(Mones.array()))
        #self.sqrtM.set_diagonal(s)
             
        self.R = _BilaplacianR(self.A, self.Msolver)      
        self.Rsolver = _BilaplacianRsolver(self.Asolver, self.M)
         
        self.mean = mean
        
        if self.mean is None:
            self.mean = dl.Vector()
            self.init_vector(self.mean, 0)
     
    def init_vector(self,x,dim):
        if dim == "noise":
            self.sqrtM.init_vector(x, 1)
        else:
            self.A.init_vector(x,dim)   
        
    def sample(self, noise, s, add_mean=True):
        rhs = self.sqrtM*noise
        self.Asolver.solve(s, rhs)
        
        if add_mean:
            s.axpy(1., self.mean)
        
class ConstrainedBiLaplacianPrior(_Prior):
    
    def __init__(self, Vh, gamma, delta, locations, m_true, anis_diff = None, pen = 1e4, rel_tol=1e-12, max_iter=1000):
        
        assert delta != 0. or pen != 0, "Intrinsic Gaussian Prior are not supported"
        self.Vh = Vh
        
        trial = dl.TrialFunction(Vh)
        test  = dl.TestFunction(Vh)
        
        if anis_diff == None:
            varfL = dl.inner(dl.nabla_grad(trial), dl.nabla_grad(test))*dl.dx
        else:
            varfL = dl.inner(anis_diff*dl.grad(trial), dl.grad(test))*dl.dx
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
        self.Asolver = dl.PETScKrylovSolver("cg", "amg")
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
        if dim == "noise":
            self.sqrtM.init_vector(x, 1)
        else:
            self.A.init_vector(x,dim)   
        
    def sample(self, noise, s, add_mean=True):
        rhs = self.sqrtM*noise
        self.Asolver.solve(s, rhs)
        
        if add_mean:
            s.axpy(1., self.mean)

class MollifiedBiLaplacianPrior(_Prior):
    
    def __init__(self, Vh, gamma, delta, locations, m_true, anis_diff = None, pen = 1e4, order=2, rel_tol=1e-12, max_iter=1000):
        
        assert delta != 0. or pen != 0, "Intrinsic Gaussian Prior are not supported"
        self.Vh = Vh
        
        trial = dl.TrialFunction(Vh)
        test  = dl.TestFunction(Vh)
        
        if anis_diff == None:
            varfL = dl.inner(dl.nabla_grad(trial), dl.nabla_grad(test))*dl.dx
        else:
            varfL = dl.inner(anis_diff*dl.grad(trial), dl.grad(test))*dl.dx
        varfM = dl.inner(trial,test)*dl.dx
        
        self.M = dl.assemble(varfM)
        self.Msolver = dl.PETScKrylovSolver("cg", "jacobi")
        self.Msolver.set_operator(self.M)
        self.Msolver.parameters["maximum_iterations"] = max_iter
        self.Msolver.parameters["relative_tolerance"] = rel_tol
        self.Msolver.parameters["error_on_nonconvergence"] = True
        self.Msolver.parameters["nonzero_initial_guess"] = False
        
        #mfun = Mollifier(gamma/delta, dl.inv(anis_diff), order, locations)
        mfun = dl.Expression(code_Mollifier)
        mfun.l = gamma/delta
        mfun.o = order
        mfun.theta0 = 1./anis_diff.theta0
        mfun.theta1 = 1./anis_diff.theta1
        mfun.alpha = anis_diff.alpha
        for ii in range(locations.shape[0]):
            mfun.addLocation(locations[ii,0], locations[ii,1])
            
        varfmo = mfun*dl.inner(trial,test)*dl.dx
        MO = dl.assemble(pen*varfmo)
                
        self.A = dl.assemble(gamma*varfL+delta*varfM + pen*varfmo)
     
        self.Asolver = dl.PETScKrylovSolver("cg", "amg")
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
        if dim == "noise":
            self.sqrtM.init_vector(x, 1)
        else:
            self.A.init_vector(x,dim)   
        
    def sample(self, noise, s, add_mean=True):
        rhs = self.sqrtM*noise
        self.Asolver.solve(s, rhs)
        
        if add_mean:
            s.axpy(1., self.mean)


        