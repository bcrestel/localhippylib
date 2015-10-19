class PDEProblem:
    """ Consider the PDE Problem:
        Given a, find u s.t. 
        F(u,a,p) = ( f(u,a), p) = 0 for all p.
        Here F is linear in p, but it may be non-linear in u and a.
    """
        
    def generate_state(self):
        """ return a vector in the shape of the state """
    
    def solveFwd(self, state, x):
        """ Solve the possibly nonlinear Fwd Problem:
        Given a, find u such that
        \delta_p F(u,a,p;\hat_p) = 0 \for all \hat_p"""
        
    def setLinearizationPoint(self,x):
        """ Set the values of the state and parameter
            for the incremental Fwd and Adj solvers """
        
    def solveIncremental(self, out, rhs, is_adj):
        """ If is_adj = False:
            Solve the forward incremental system:
            Given u, a, find \tilde_u s.t.:
            \delta_{pu} F(u,a,p; \hat_p, \tilde_u) = rhs for all \hat_p.
            
            If is_adj = True:
            Solve the adj incremental system:
            Given u, a, find \tilde_p s.t.:
            \delta_{up} F(u,a,p; \hat_u, \tilde_p) = rhs for all \delta_u.
        """
        
    def applyF_pa(self, dir, out, is_adj):
        """ If is_adj = False:
            Given u, a; compute 
            \delta_{pa} F(u,a,p; \hat_p, \tilde_a) in the direction \tilde_a = dir for all \hat_p

            If is_adj = True:
            Given u, a; compute \delta_{ap} F(u,a,p; \hat_a, \tilde_p) in the direction \tilde_p = dir for all \hat_a
        """
        
    def applyF_aa(self, delta_a, out):
        """ 
        Given u,a,p; compute \delta_{aa} F(u,a,p; \hat_a, \tilde_a) in the direction \tilde_a = dir for all \hat_a
        """
        
    def applyF_au(self, dir, out, is_adj):
        """ If is_adj = False:
            Given u, a, p; compute \delta_{au} F(u,a,p; \hat_a, tilde_u) in the direction \tilde_u = dir for all \hat_a
            
            If is_adj = True:
            Given u, a, p; compute \delta_{ua} F(u,a,p; \hat_u, tilde_a) in the direction \tilde_a = dir for all \hat_u
        """

import dolfin as dl
import math
import sys
sys.path.append( "../../" )
from hippylib import *

class IceProblem(PDEProblem):
    
    def __init__(self, Vh, ds_robin, RateFactor, GlenExp, density, gravity, eps, pen=1e3):
        """ Construct an IceProblem with the following paramenters:
        - Vh: Finite Element Space for the state variable
        - ds_robin: Boundary integrate for the Robin Boundary Condition
        - RateFactor: The (temperature dependent) flow rate factor 
        - GlenExp: The Glen's flow law exponent
        - density: The ice density
        - gravity: The gravity field.
        - eps: a regularization parameter for the Glen Flow effective viscosity. 
        """
        
        self.Vh = Vh
        self.ds_robin = ds_robin
        self.RateFactor = RateFactor
        self.GlenExp = GlenExp
        self.density = density
        self.gravity = gravity
        self.eps = eps
        self.pen = pen
        
        up = dl.TrialFunction(self.Vh[STATE])
        vq = dl.TestFunction(self.Vh[STATE])
        self.M = dl.assemble( dl.inner(up,vq)*dl.dx )
        self.Msolver = dl.PETScKrylovSolver("cg", "jacobi")
        self.Msolver.set_operator(self.M)
        
        self.newton_parameters = {}
        self.newton_parameters["rel_tolerance"] = 1e-6
        self.newton_parameters["abs_tolerance"] = 1e-9
        self.newton_parameters["gdu_tolerance"] = 1e-18
        self.newton_parameters["lin_tolerance"] = 1e-12
        self.newton_parameters["max_iter"]      = 100
        self.newton_parameters["c_armijo"]      = 1e-5
        self.newton_parameters["max_backtracking_iter"] = 10
        
    def generate_state(self):
        return dl.Function(self.Vh[STATE]).vector()
    
    def normF(self, F):
        f = dl.Vector()
        self.M.init_vector(f,1)
        self.Msolver.solve(f,F)
        return math.sqrt( f.inner(F) )
    
    def FNewtonean(self,up,a,vq, eta_noA):
        n = dl.FacetNormal(self.Vh[STATE].mesh())
        h = dl.CellSize(self.Vh[STATE].mesh())
        
        uh,ph = dl.split(up)
        vh,qh = dl.split(vq)
        
        myA = self.RateFactor**(-1./self.GlenExp)
        eta = myA*eta_noA
                
        def t(u,p,n): return dl.dot(eta*dl.sym(dl.nabla_grad(u)),n) - p*n
        
        a11 = eta*dl.inner(dl.sym( dl.nabla_grad(uh) ), dl.sym( dl.nabla_grad(vh) ))*dl.dx
        a12 = -dl.nabla_div(vh)*ph*dl.dx
        a21 = -dl.nabla_div(uh)*qh*dl.dx
        robinvarf  = dl.exp(a)*dl.inner(uh - dl.dot(dl.outer(n,n),uh), vh) * self.ds_robin
        weak_bc = -dl.dot(n, t(uh, ph, n) )*dl.dot(vh,n)*self.ds_robin \
                  - dl.dot(n, t(vh, qh, n) )*dl.dot(uh,n)*self.ds_robin \
                  + self.pen/h*dl.dot(uh,n)*dl.dot(vh, n)*self.ds_robin
        #weak_bc = pen/h*dl.dot(uh,n)*dl.dot(vh, n)*self.ds_robin
        
        f1 = dl.inner(vh, self.gravity*self.density)*dl.dx
        f2 = dl.Constant(0.)*qh*dl.dx

        if type(up) == dl.Function:        
            return a11+a12+a21+robinvarf+weak_bc - f1-f2
        elif type(up) == dl.Argument:
            return a11+a12+a21+robinvarf+weak_bc
        else:
            raise Exception()

    def Energy(self,up, a):

        n = dl.FacetNormal(self.Vh[STATE].mesh())
        h = dl.CellSize(self.Vh[STATE].mesh())
        
        uh,ph = dl.split(up)
        
        myA = self.RateFactor**(-1./self.GlenExp)
        myExp = dl.Constant(.5)*( dl.Constant(1.) + self.GlenExp)/self.GlenExp
        
        Phi = dl.Constant(1.)/myExp*myA*( dl.Constant(.5)*dl.inner( dl.sym(dl.nabla_grad(uh)), dl.sym(dl.nabla_grad(uh)) )  + self.eps )**myExp
                
        a11 = Phi*dl.dx
        robinvarf  = dl.Constant(.5)*dl.exp(a)*dl.inner(uh - dl.dot(dl.outer(n,n),uh), uh) * self.ds_robin
        weak_bc = dl.Constant(.5)*self.pen/h*dl.dot(uh,n)*dl.dot(uh, n)*self.ds_robin
        
        f1 = dl.inner(uh, self.gravity*self.density)*dl.dx
       
        return a11+robinvarf+weak_bc - f1
        

                
    def F(self,up,a,vq):
        n = dl.FacetNormal(self.Vh[STATE].mesh())
        h = dl.CellSize(self.Vh[STATE].mesh())
        
        uh,ph = dl.split(up)
        vh,qh = dl.split(vq)
        
        myA = self.RateFactor**(-1./self.GlenExp)
        myExp = dl.Constant(.5)*( dl.Constant(1.) - self.GlenExp)/self.GlenExp
        
        eta = myA*( dl.Constant(.5)*dl.inner( dl.sym(dl.nabla_grad(uh)), dl.sym(dl.nabla_grad(uh)) )  + self.eps )**myExp
        #eta = myA
        
        def t(u,p,n): return dl.dot(eta*dl.sym(dl.nabla_grad(u)),n) - p*n
        
        a11 = eta*dl.inner(dl.sym( dl.nabla_grad(uh) ), dl.sym( dl.nabla_grad(vh) ))*dl.dx
        a12 = -dl.nabla_div(vh)*ph*dl.dx
        a21 = -dl.nabla_div(uh)*qh*dl.dx
        robinvarf  = dl.exp(a)*dl.inner(uh - dl.dot(dl.outer(n,n),uh), vh) * self.ds_robin
        #weak_bc = -dl.dot(n, t(uh, ph, n) )*dl.dot(vh,n)*self.ds_robin \
        #          - dl.dot(n, t(vh, qh, n) )*dl.dot(uh,n)*self.ds_robin \
        #          + self.pen/h*dl.dot(uh,n)*dl.dot(vh, n)*self.ds_robin
        weak_bc = self.pen/h*dl.dot(uh,n)*dl.dot(vh, n)*self.ds_robin
        
        f1 = dl.inner(vh, self.gravity*self.density)*dl.dx
        f2 = dl.Constant(0.)*qh*dl.dx

        if type(up) == dl.Function:        
            return a11+a12+a21+robinvarf+weak_bc - f1-f2
        elif type(up) == dl.Argument:
            return a11+a12+a21+robinvarf+weak_bc
        else:
            raise Exception()
        
    def solveFwd(self, state, x):
        # Assume u0 is div free
        up = dl.Function(self.Vh[STATE], state)
        a  = dl.Function(self.Vh[PARAMETER], x[PARAMETER])
        vq = dl.TestFunction(self.Vh[STATE])
        
        E0 = dl.assemble( self.Energy(up, a) )
        F0 = dl.assemble( self.F(up, a, vq) )
        norm_F0 = self.normF(F0)
        
        rtol = self.newton_parameters["rel_tolerance"]
        atol = self.newton_parameters["abs_tolerance"]
        gdu_tol  = self.newton_parameters["gdu_tolerance"]
        lin_tol = self.newton_parameters["lin_tolerance"]
        maxiter = self.newton_parameters["max_iter"]
        c_armijo = self.newton_parameters["c_armijo"]
        maxiter_backtracking = self.newton_parameters["max_backtracking_iter"]
        converged = False
        
        tol = max( norm_F0*rtol, atol)
        
        up_prev = dl.Function(self.Vh[STATE])
        du = self.generate_state()
        for it in range(maxiter):
            f_form = self.F(up,a,vq)
            J_form = dl.derivative(f_form, up)
            J, F = dl.assemble_system(J_form, f_form)
            
            norm_F = self.normF(F)
            
            if norm_F < tol:
                converged = True
                break
            
            Jsolver = dl.PETScLUSolver(J)
            Jsolver.solve(du, F)
            
            duF = F.inner(du)
            
            if duF < gdu_tol:
                converged = True
                break
            
            bt_converged = False
            alpha = 1.
            E_prev = dl.assemble( self.Energy(up, a) )
            up_prev.assign(up)
            
            for j in range(maxiter_backtracking):
                up.vector().axpy(-alpha, du)
                E_cur = dl.assemble( self.Energy(up, a) )
                if E_cur < E_prev - alpha*c_armijo*duF:
                    bt_converged = True
                    break
                else:
                    up.assign(up_prev)
                    alpha *= .5
                    
            if not bt_converged:
                break
                        
        assert converged
                    
                
            
            
