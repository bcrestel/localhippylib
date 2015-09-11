import dolfin as dl
import sys
sys.path.append( "../../" )
from hippylib import *
import numpy as np
import matplotlib.pyplot as plt

def u_boundary(x, on_boundary):
    return on_boundary and ( x[1] < dl.DOLFIN_EPS or x[1] > 1.0 - dl.DOLFIN_EPS)

def v_boundary(x, on_boundary):
    return on_boundary and ( x[0] < dl.DOLFIN_EPS or x[0] > 1.0 - dl.DOLFIN_EPS)

def compute_velocity(mesh, Vh, a, u):
    #export the velocity field v = - exp( a ) \grad u: then we solve ( exp(-a) v, w) = ( u, div w)
    Vv = dl.FunctionSpace(mesh, 'RT', 1)
    v = dl.Function(Vv, name="velocity")
    vtrial = dl.TrialFunction(Vv)
    vtest = dl.TestFunction(Vv)
    afun = dl.Function(Vh[CONTROL], a)
    ufun = dl.Function(Vh[STATE], u)
    Mv = dl.exp(-afun) *dl.inner(vtrial, vtest) *dl.dx
    n = dl.FacetNormal(mesh)
    class TopBoundary(dl.SubDomain):
        def inside(self,x,on_boundary): 
            return on_boundary and x[1] > 1 - dl.DOLFIN_EPS
        
    Gamma_M = TopBoundary()
    boundaries = dl.FacetFunction("size_t", mesh)
    boundaries.set_all(0)
    Gamma_M.mark(boundaries, 1)
    dss = dl.Measure("ds")[boundaries]
    rhs = ufun*dl.div(vtest)*dl.dx - dl.dot(vtest,n)*dss(1)
    bcv = dl.DirichletBC(Vv, dl.Expression( ("0.0", "0.0") ), v_boundary)
    dl.solve(Mv == rhs, v, bcv)
    
    return v

def true_model(Vh, gamma, delta, anis_diff):
    prior = BiLaplacianPrior(Vh, gamma, delta, anis_diff )
    noise = dl.Vector()
    prior.init_vector(noise,"noise")
    noise_size = noise.array().shape[0]
    noise.set_local( np.random.randn( noise_size ) )
    atrue = dl.Vector()
    prior.init_vector(atrue, 0)
    prior.sample(noise,atrue)
    return atrue

class Poisson:
    def __init__(self, mesh, Vh, targets, prior):
        """
        Construct a model by proving
        - the mesh
        - the finite element spaces for the STATE/ADJOINT variable and the control variable
        - the Prior information
        """
        self.mesh = mesh
        self.Vh = Vh
        
        # Initialize Expressions
        self.f = dl.Expression("0.0")
        
        self.u_bdr = dl.Expression("x[1]")
        self.u_bdr0 = dl.Expression("0.0")
        self.bc = dl.DirichletBC(self.Vh[STATE], self.u_bdr, u_boundary)
        self.bc0 = dl.DirichletBC(self.Vh[STATE], self.u_bdr0, u_boundary)
                
        # Assemble constant matrices      
        self.M = self.assembleM()
        self.prior = prior
        self.B = assemblePointwiseObservation(Vh[STATE],targets)
                
        self.A = []
        self.At = []
        self.C = []
        self.Raa = []
        self.Wau = []
        
        self.u_o = dl.Vector()
        self.B.init_vector(self.u_o,0)
        self.noise_variance = 0
        
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
            x = [dl.Vector(), dl.Vector(), dl.Vector()]
            self.B.init_vector(x[STATE],1)
            self.prior.init_vector(x[CONTROL],0)
            self.B.init_vector(x[ADJOINT], 1)
        elif component == STATE:
            x = dl.Vector()
            self.B.init_vector(x,1)
        elif component == CONTROL:
            x = dl.Vector()
            self.prior.init_vector(x,0)
        elif component == ADJOINT:
            x = dl.Vector()
            self.B.init_vector(x,1)
            
        return x
    
    def init_control(self, a):
        """
        Reshape a so that it is compatible with the control variable
        """
        self.prior.init_vector(a,0)
        
    def assembleA(self,x, assemble_adjoint = False, assemble_rhs = False):
        """
        Assemble the matrices and rhs for the forward/adjoint problems
        """
        trial = dl.TrialFunction(self.Vh[STATE])
        test = dl.TestFunction(self.Vh[STATE])
        c = dl.Function(self.Vh[CONTROL], x[CONTROL])
        Avarf = dl.inner(dl.exp(c)*dl.nabla_grad(trial), dl.nabla_grad(test))*dl.dx
        if not assemble_adjoint:
            bform = dl.inner(self.f, test)*dl.dx
            Matrix, rhs = dl.assemble_system(Avarf, bform, self.bc)
        else:
            # Assemble the adjoint of A (i.e. the transpose of A)
            s = dl.Function(self.Vh[STATE], x[STATE])
            bform = dl.inner(dl.Constant(0.), test)*dl.dx
            Matrix, _ = dl.assemble_system(dl.adjoint(Avarf), bform, self.bc0)
            Bu = -(self.B*x[STATE])
            Bu += self.u_o
            rhs = dl.Vector()
            self.B.init_vector(rhs, 1)
            self.B.transpmult(Bu,rhs)
            rhs *= 1.0/self.noise_variance
            
        if assemble_rhs:
            return Matrix, rhs
        else:
            return Matrix
    
    def assembleC(self, x):
        """
        Assemble the derivative of the forward problem with respect to the control
        """
        trial = dl.TrialFunction(self.Vh[CONTROL])
        test = dl.TestFunction(self.Vh[STATE])
        s = dl.Function(Vh[STATE], x[STATE])
        c = dl.Function(Vh[CONTROL], x[CONTROL])
        Cvarf = dl.inner(dl.exp(c) * trial * dl.nabla_grad(s), dl.nabla_grad(test)) * dl.dx
        C = dl.assemble(Cvarf)
#        print "||c||", x[CONTROL].norm("l2"), "||s||", x[STATE].norm("l2"), "||C||", C.norm("linf")
        self.bc0.zero(C)
        return C
        
    def assembleM(self):
        """
        Assemble the mass matrix in the control space.
        This is needed in evalGradientControl to compute the L2 norm of the gradient
        """
        trial = dl.TrialFunction(self.Vh[CONTROL])
        test  = dl.TestFunction(self.Vh[CONTROL])
        varf = dl.inner(trial, test)*dl.dx
        return dl.assemble(varf)
        
    def assembleWau(self, x):
        """
        Assemble the derivative of the control equation with respect to the state
        """
        trial = dl.TrialFunction(self.Vh[STATE])
        test  = dl.TestFunction(self.Vh[CONTROL])
        a = dl.Function(self.Vh[ADJOINT], x[ADJOINT])
        c = dl.Function(self.Vh[CONTROL], x[CONTROL])
        varf = dl.inner(dl.exp(c)*dl.nabla_grad(trial),dl.nabla_grad(a))*test*dl.dx
        Wau = dl.assemble(varf)
        dummy = dl.Vector()
        Wau.init_vector(dummy,0)
        self.bc0.zero_columns(Wau, dummy)
        return Wau
    
    def assembleRaa(self, x):
        """
        Assemble the derivative of the control equation with respect to the control (Newton method)
        """
        trial = dl.TrialFunction(self.Vh[CONTROL])
        test  = dl.TestFunction(self.Vh[CONTROL])
        s = dl.Function(self.Vh[STATE], x[STATE])
        c = dl.Function(self.Vh[CONTROL], x[CONTROL])
        a = dl.Function(self.Vh[ADJOINT], x[ADJOINT])
        varf = dl.inner(dl.nabla_grad(a),dl.exp(c)*dl.nabla_grad(s))*trial*test*dl.dx
        return dl.assemble(varf)

            
    def cost(self, x):
        """
        Given the list x = [u,a,p] which describes the state, control, and
        adjoint variable compute the cost functional as the sum of 
        the misfit functional and the regularization functional.
        
        Return the list [cost functional, regularization functional, misfit functional]
        
        Note: p is not needed to compute the cost functional
        """        
        assert x[STATE] != None
                       
        diff = self.B*x[STATE]
        diff -= self.u_o
        misfit = (.5/self.noise_variance) * diff.inner(diff)
        
        Rdiff_x = dl.Vector()
        self.prior.init_vector(Rdiff_x,0)
        diff_x = x[CONTROL] - self.prior.mean
        self.prior.R.mult(diff_x, Rdiff_x)
        reg = .5 * diff_x.inner(Rdiff_x)
        
        c = misfit + reg
        
        return c, reg, misfit
    
    def solveFwd(self, out, x, tol=1e-9):
        """
        Solve the forward problem.
        """
        A, b = self.assembleA(x, assemble_rhs = True)
        A.init_vector(out, 1)
        solver = dl.PETScKrylovSolver("cg", amg_method())
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
        
        solver = dl.PETScKrylovSolver("cg", amg_method())
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

        self.prior.init_vector(mg,0)
        C.transpmult(x[ADJOINT], mg)
        Rdx = dl.Vector()
        self.prior.init_vector(Rdx,0)
        dx = x[CONTROL] - self.prior.mean
        self.prior.R.mult(dx, Rdx)   
        mg.axpy(1., Rdx)
        
        g = dl.Vector()
        self.prior.init_vector(g,1)
        
        solver = dl.PETScKrylovSolver("cg", "jacobi")
        solver.set_operator(self.M)
        solver.solve(g, mg)
        g_norm = dl.sqrt( g.inner(mg) )
        
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
        solver = dl.PETScKrylovSolver("cg", amg_method())
        solver.set_operator(self.A)
        solver.parameters["relative_tolerance"] = tol
        self.A.init_vector(sol,1)
        nit = solver.solve(sol,rhs)
#        print "FwdInc", (self.A*sol-rhs).norm("l2")/rhs.norm("l2"), nit
        
    def solveAdjIncremental(self, sol, rhs, tol):
        """
        Solve the incremental adjoint problem for a given rhs
        """
        solver = dl.PETScKrylovSolver("cg", amg_method())
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
        help = dl.Vector()
        self.B.init_vector(help, 0)
        self.B.mult(du, help)
        self.B.transpmult(help, out)
        out *= 1./self.noise_variance
    
    def applyWua(self, da, out):
        self.Wau.transpmult(da,out)

    
    def applyWau(self, du, out):
        self.Wau.mult(du, out)
    
    def applyR(self, da, out):
        self.prior.R.mult(da, out)
        
    def Rsolver(self):        
        return self.prior.Rsolver
    
    def applyRaa(self, da, out):
        self.Raa.mult(da, out)
            
if __name__ == "__main__":
    dl.set_log_active(False)
    sep = "\n"+"#"*80+"\n"
    print sep, "Set up the mesh and finite element spaces", sep
    ndim = 2
    nx = 64
    ny = 64
    mesh = dl.UnitSquareMesh(nx, ny)
    Vh2 = dl.FunctionSpace(mesh, 'Lagrange', 2)
    Vh1 = dl.FunctionSpace(mesh, 'Lagrange', 1)
    Vh = [Vh2, Vh1, Vh2]
    print "Number of dofs: STATE={0}, CONTROL={1}, ADJOINT={2}".format(Vh[STATE].dim(), Vh[CONTROL].dim(), Vh[ADJOINT].dim())
    
    print sep, "Set up the location of observation, Prior Information, and model", sep
    ntargets = 300
    np.random.seed(seed=1)
    targets = np.random.uniform(0.1,0.9, [ntargets, ndim] )
    print "Number of observation points: {0}".format(ntargets)
    
    gamma = .1
    delta = .5
    
#    myangle = dl.Constant(math.pi/4)
#    mydir1 = dl.Constant(2.)
#    mydir2 = dl.Constant(.5)
#    c00 = mydir1*dl.sin(myangle)*dl.sin(myangle) + mydir2*dl.cos(myangle)*dl.cos(myangle)
#    c11 = mydir2*dl.sin(myangle)*dl.sin(myangle) + mydir1*dl.cos(myangle)*dl.cos(myangle)
#    c01 = (mydir1-mydir2)*dl.sin(myangle)*dl.cos(myangle)
#    anis_diff = dl.as_matrix(((c00, c01), (c01, c11)))
    anis_diff = dl.Expression(code_AnisTensor2D)
    anis_diff.theta0 = 2.
    anis_diff.theta1 = .5
    anis_diff.alpha = math.pi/4
    atrue = true_model(Vh[CONTROL], gamma, delta,anis_diff)
        
    locations = np.array([[0.1, 0.1], [0.1, 0.9], [.5,.5], [.9, .1], [.9, .9]])
    if 1:
        pen = 1e1
        prior = MollifiedBiLaplacianPrior(Vh[CONTROL], gamma, delta, locations, atrue, anis_diff, pen)
    else:
        pen = 1e4
        prior = ConstrainedBiLaplacianPrior(Vh[CONTROL], gamma, delta, locations, atrue, anis_diff, pen)
        
    print "Prior regularization: (delta_x - gamma*Laplacian)^order: delta={0}, gamma={1}, order={2}".format(delta, gamma,2)    
            
    model = Poisson(mesh, Vh, targets, prior)
    
    #Generate synthetic observations
    x = [model.generate_vector(STATE), atrue, None]
    model.solveFwd(x[STATE], x, 1e-9)
    utrue = x[STATE]
    model.B.mult(x[STATE], model.u_o)
    rel_noise = 0.01
    MAX = model.u_o.norm("linf")
    noise_std_dev = rel_noise * MAX
    randn_perturb(model.u_o, noise_std_dev)
    model.noise_variance = noise_std_dev*noise_std_dev
   
    print sep, "Test the gradient and the Hessian of the model", sep
    a0 = dl.interpolate(dl.Expression("sin(x[0])"), Vh[CONTROL])
    modelVerify(model, a0.vector(), 1e-4, 1e-6)

    print sep, "Find the MAP point", sep
    a0 = prior.mean.copy()
    solver = ReducedSpaceNewtonCG(model)
    solver.parameters["rel_tolerance"] = 1e-9
    solver.parameters["abs_tolerance"] = 1e-12
    solver.parameters["max_iter"]      = 25
    solver.parameters["inner_rel_tolerance"] = 1e-15
    solver.parameters["c_armijo"] = 1e-4
    solver.parameters["GN_iter"] = 5
    
    x = solver.solve(a0)
    
    if solver.converged:
        print "\nConverged in ", solver.it, " iterations."
    else:
        print "\nNot Converged"

    print "Termination reason: ", solver.termination_reasons[solver.reason]
    print "Final gradient norm: ", solver.final_grad_norm
    print "Final cost: ", solver.final_cost
        
    print sep, "Compute the low rank Gaussian Approximation of the posterior", sep
    model.setPointForHessianEvaluations(x)
    Hmisfit = ReducedHessian(model, solver.parameters["inner_rel_tolerance"], gauss_newton_approx=False, misfit_only=True)
    k = 50
    p = 20
    print "Double Pass Algorithm. Requested eigenvectors: {0}; Oversampling {1}.".format(k,p)
    Omega = np.random.randn(x[CONTROL].array().shape[0], k+p)
    #d, U = singlePassG(Hmisfit, model.R, model.Rsolver, Omega, k, check_Bortho=True, check_Aortho=True, check_residual=True)
    d, U = doublePassG(Hmisfit, prior.R, prior.Rsolver, Omega, k, check_Bortho=False, check_Aortho=False, check_residual=False)
    posterior = GaussianLRPosterior(prior, d, U)
    posterior.mean = x[CONTROL]
    
    post_tr, prior_tr, corr_tr = posterior.trace(method="Estimator", tol=1e-1, min_iter=20, max_iter=100)
    print "Posterior trace {0:5e}; Prior trace {1:5e}; Correction trace {2:5e}".format(post_tr, prior_tr, corr_tr)
    post_pw_variance, pr_pw_variance, corr_pw_variance = posterior.pointwise_variance("Exact")

    print sep, "Save State, Control, Adjoint, and observation in paraview", sep
    xxname = ["State", "Control", "Adjoint"]
    xx = [dl.Function(Vh[i], x[i], name=xxname[i]) for i in range(len(Vh))]
    dl.File("results/poisson_state.pvd") << xx[STATE]
    dl.File("results/poisson_state_true.pvd") << dl.Function(Vh[STATE], utrue, name = xxname[STATE])
    dl.File("results/poisson_control.pvd") << xx[CONTROL]
    dl.File("results/poisson_control_true.pvd") << dl.Function(Vh[CONTROL], atrue, name = xxname[CONTROL])
    dl.File("results/poisson_control_prmean.pvd") << dl.Function(Vh[CONTROL], prior.mean, name = xxname[CONTROL])
    dl.File("results/poisson_adjoint.pvd") << xx[ADJOINT]
    
    vtrue = compute_velocity(mesh, Vh, atrue, utrue)
    dl.File("results/poisson_vel_true.pvd") << vtrue
    v_map = compute_velocity(mesh, Vh, x[CONTROL], x[STATE])
    dl.File("results/poisson_vel.pvd") << v_map
    
    exportPointwiseObservation(targets, model.u_o, "results/poisson_observation.vtp")
    
    fid = dl.File("results/pointwise_variance.pvd")
    fid << dl.Function(Vh[CONTROL], post_pw_variance, name="Posterior")
    fid << dl.Function(Vh[CONTROL], pr_pw_variance, name="Prior")
    fid << dl.Function(Vh[CONTROL], corr_pw_variance, name="Correction")
    
    
    print sep, "Generate samples from Prior and Posterior\n","Export generalized Eigenpairs", sep
    fid_prior = dl.File("samples/sample_prior.pvd")
    fid_post  = dl.File("samples/sample_post.pvd")
    nsamples = 500
    noise = dl.Vector()
    posterior.init_vector(noise,"noise")
    noise_size = noise.array().shape[0]
    s_prior = dl.Function(Vh[CONTROL], name="sample_prior")
    s_post = dl.Function(Vh[CONTROL], name="sample_post")
    for i in range(nsamples):
        noise.set_local( np.random.randn( noise_size ) )
        posterior.sample(noise, s_prior.vector(), s_post.vector())
        fid_prior << s_prior
        fid_post << s_post
        
    #Save eigenvalues for printing:
    posterior.exportU(Vh[CONTROL], "hmisfit/evect.pvd")
    np.savetxt("hmisfit/eigevalues.dat", d)
    
    print sep, "Visualize results", sep
    dl.plot(xx[STATE], title = xxname[STATE])
    dl.plot(dl.exp(xx[CONTROL]), title = xxname[CONTROL])
    dl.plot(xx[ADJOINT], title = xxname[ADJOINT])
    
    plt.plot(range(0,k), d, 'b*', range(0,k), np.ones(k), '-r')
    plt.yscale('log')
        
    plt.show()    
    dl.interactive()
    
