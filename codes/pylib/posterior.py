from dolfin import Vector, Function, File
from lowRankOperator import LowRankOperator
import numpy as np

class LowRankHessian:
    def __init__(self, prior, d, U):
        self.prior = prior
        self.LowRankH = LowRankOperator(d, U)
        dsolve = d / (np.ones(d.shape, dtype=d.dtype) + d)
        self.LowRankHinv = LowRankOperator(dsolve, U)
        self.help = Vector()
        self.init_vector(self.help, 0)
        
    def init_vector(self,x, dim):
        self.prior.init_vector(x,dim)
        
    def mult(self, x, y):
        self.prior.R.mult(x,y)
        self.LowRankH.mult(x, self.help)
        y.axpy(1, self.help)
        
        
    def solve(self, sol, rhs):
        self.prior.Rsolver.solve(sol, rhs)
        self.LowRankHinv.mult(rhs, self.help)
        sol.axpy(-1, self.help)
        
class LowRankPosteriorSampler:
    def __init__(self, prior, d, U):
        self.prior = prior        
        ones = np.ones( d.shape, dtype=d.dtype )        
        self.d = np.power(ones + d, -.5) - ones
        self.lrsqrt = LowRankOperator(self.d, U)
        self.help = Vector()
        self.init_vector(self.help, 0)
        
    def init_vector(self,x, dim):
        self.prior.init_vector(x,dim)
        
    def sample(self, noise, s):
        self.prior.R.mult(noise, self.help)
        self.lrsqrt.mult(self.help, s)
        s.axpy(1, noise)

class GaussianLRPosterior:
    def __init__(self, prior, d, U, mean=None):
        self.prior = prior
        self.d = d
        self.U = U
        self.Hlr = LowRankHessian(prior, d, U)
        self.sampler = LowRankPosteriorSampler(self.prior, self.d, self.U)
        self.mean=None
        
    def init_vector(self,x, dim):
        self.prior.init_vector(x,dim)
        
    def sample(self, *args, **kwargs):
        """
        possible calls:
        1) sample(s_prior, s_post, add_mean=True)
           - s_prior is a sample from the prior centered at 0 (input)
           - s_post is a sample from the posterior (output)
           - if add_mean=True (default) than the samples will be centered at the map
             point
        2) sample(noise, s_prior, s_post, add_mean=True)
           - noise is a realization of white noise (input)
           - s_prior is a sample from the prior (output)
           - s_post  is a sample from the posterior
           - if add_mean=True (default) than the prior and posterior samples will be
                centered at the respective means.
        """
        add_mean = True
        for name, value in kwargs.items():
            if name == "add_mean":
                add_mean = value
            else:
                raise NameError(name)
        
        if len(args) == 2:
            self._sample_given_prior(args[0], args[1])
            if add_mean:
                args[1].axpy(1., self.mean)
        elif len(args) == 3:
            self._sample_given_white_noise(args[0], args[1], args[2])
            if add_mean:
                    args[1].axpy(1., self.prior.mean) 
                    args[2].axpy(1., self.mean)
        else:
            raise NameError('Invalid number of parameters in Posterior::sample')
        
    def _sample_given_white_noise(self, noise, s_prior, s_post):
        self.prior.sample(noise, s_prior, add_mean=False)
        self.sampler.sample(s_prior, s_post)
        
    def _sample_given_prior(self,s_prior, s_post):
        self.sampler.sample(s_prior, s_post)
    
    def exportU(self, Vh, fname, varname = "evect", normalize=1):
        #Save the scaled eigenvectors (normalized such that || v ||_inf = 1)
        evect = Function(Vh, name=varname)
        fid = File(fname)
        
        for i in range(0,self.U.shape[1]):
            Ui = self.U[:,i]
            if normalize:
                s = 1/np.linalg.norm(Ui, np.inf)
                evect.vector().set_local(s*Ui)
            else:
                evect.vector().set_local(Ui)
            fid << evect
            
    def trace(self, method="Exact", tol=1e-1, min_iter=20, max_iter=100):
        pr_trace = self.prior.trace(method, tol, min_iter, max_iter)
        corr_trace = self.Hlr.LowRankHinv.trace(self.prior.M)
        post_trace = pr_trace - corr_trace
        return post_trace, pr_trace, corr_trace
    
    def pointwise_variance(self, method="Exact", path_len = 12):
        pr_pointwise_variance = self.prior.pointwise_variance(method, path_len)
        correction_pointwise_variance = Vector()
        self.init_vector(correction_pointwise_variance, 0)
        self.Hlr.LowRankHinv.get_diagonal(correction_pointwise_variance)
        post_pointwise_variance = pr_pointwise_variance - correction_pointwise_variance
        return post_pointwise_variance, pr_pointwise_variance, correction_pointwise_variance
        
        
