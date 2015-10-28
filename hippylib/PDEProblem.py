'''
Created on Oct 19, 2015

@author: uvilla
'''
class PDEProblem:
    """ Consider the PDE Problem:
        Given a, find u s.t. 
        F(u,a,p) = ( f(u,a), p) = 0 for all p.
        Here F is linear in p, but it may be non-linear in u and a.
    """
        
    def generate_state(self):
        """ return a vector in the shape of the state """
    
    def solveFwd(self, state, x, tol):
        """ Solve the possibly nonlinear Fwd Problem:
        Given a, find u such that
        \delta_p F(u,a,p;\hat_p) = 0 \for all \hat_p"""
        
    def solveAdj(self, state, x, adj_rhs, tol):
        """ Solve the linear Adj Problem: 
            Given a, u; find p such that
            \delta_u F(u,a,p;\hat_u) = 0 \for all \hat_u
        """
     
    def eval_da(self, x, out):
        """Given u,a,p; eval \delta_a F(u,a,p; \hat_a) \for all \hat_a """
         
    def setLinearizationPoint(self,x):
        """ Set the values of the state and parameter
            for the incremental Fwd and Adj solvers """
        
    def solveIncremental(self, out, rhs, is_adj, mytol):
        """ If is_adj = False:
            Solve the forward incremental system:
            Given u, a, find \tilde_u s.t.:
            \delta_{pu} F(u,a,p; \hat_p, \tilde_u) = rhs for all \hat_p.
            
            If is_adj = True:
            Solve the adj incremental system:
            Given u, a, find \tilde_p s.t.:
            \delta_{up} F(u,a,p; \hat_u, \tilde_p) = rhs for all \delta_u.
        """
    
    def apply_ij(self,i,j, dir, out):   
        """
            Given u, a, p; compute 
            \delta_{ij} F(u,a,p; \hat_i, \tilde_j) in the direction \tilde_j = dir for all \hat_i
        """

class PDEVariationalProblem(PDEProblem):
    def __init__(self, Vh, varf_handler, bc, bc0):
        self.Vh = Vh
        self.varf_handler = varf_handler
        self.bc = bc
        self.bc0 = bc0
        
    def generate_state(self):
        """ return a vector in the shape of the state """
        return dl.Function(self.Vh[STATE]).vector()
    
    def solveFwd(self, state, x, tol):
        """ Solve the possibly nonlinear Fwd Problem:
        Given a, find u such that
        \delta_p F(u,a,p;\hat_p) = 0 \for all \hat_p"""
        state.set_local(x[STATE].array())
        u = dl.Function(self.Vh[STATE], state)
        a = dl.Function(self.Vh[PARAMETER], x[PARAMETER])
        p = dl.TestFunction(self.Vh[ADJOINT])
        res_form = self.varf_handler(u,a,p)
        dl.solve(res_form == 0, u, self.bc)
        
    def solveAdj(self, adj, x, adj_rhs, tol):
        """ Solve the linear Adj Problem: 
            Given a, u; find p such that
            \delta_u F(u,a,p;\hat_u) = 0 \for all \hat_u
        """
        u = dl.Function(self.Vh[STATE], x[STATE])
        a = dl.Function(self.Vh[PARAMETER], x[PARAMETER])
        p = dl.TrialFunction(self.Vh[ADJOINT])
        res_form = self.varf_handler(u,a,p)
        adj_form = dl.derivative(res_form, u)
        Aadj = dl.assemle(adj_form)
        self.bc0.apply(Aadj, adj_rhs)
        dl.solve(Aadj, adj, adj_rhs)
     
    def eval_da(self, x, out):
        """Given u,a,p; eval \delta_a F(u,a,p; \hat_a) \for all \hat_a """
         
    def setLinearizationPoint(self,x):
        """ Set the values of the state and parameter
            for the incremental Fwd and Adj solvers """
        
    def solveIncremental(self, out, rhs, is_adj, mytol):
        """ If is_adj = False:
            Solve the forward incremental system:
            Given u, a, find \tilde_u s.t.:
            \delta_{pu} F(u,a,p; \hat_p, \tilde_u) = rhs for all \hat_p.
            
            If is_adj = True:
            Solve the adj incremental system:
            Given u, a, find \tilde_p s.t.:
            \delta_{up} F(u,a,p; \hat_u, \tilde_p) = rhs for all \delta_u.
        """
    
    def apply_ij(self,i,j, dir, out):   
        """
            Given u, a, p; compute 
            \delta_{ij} F(u,a,p; \hat_i, \tilde_j) in the direction \tilde_j = dir for all \hat_i
        """

    