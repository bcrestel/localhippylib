import numpy as np
import dolfin as dl

from hippylib.linalg import vector2Function
from hippylib.variables import STATE, ADJOINT, PARAMETER


class JointModel:
    """ This class builds a model for joint inversion in the case of two
    independent physics connected only via the regularization term (or prior) """

    def __init__(self, model1, model2, jointregularization, alphareg=1.0):
        """ Provides models model1 and model2 with ZeroPrior as their
        regularization, and provide jointregularization. """

        self.model1 = model1
        self.model2 = model2

        # individual models shall not have their own regularization!
        assert self.model1.Prior.isZeroPrior()
        assert self.model2.Prior.isZeroPrior()

        self.Vh = [None, None, None]
        self.Vh[STATE] = self.model1.problem.Vh[STATE] * self.model2.problem.Vh[STATE]
        self.Vh[ADJOINT] = self.model1.problem.Vh[ADJOINT] * self.model2.problem.Vh[ADJOINT]
        self.Vh[PARAMETER] = self.model1.problem.Vh[PARAMETER] * self.model2.problem.Vh[PARAMETER]

        self.M = [None, None, None]
        test, trial = dl.TestFunction(self.Vh[STATE]), dl.TrialFunction(self.Vh[STATE])
        self.M[STATE] = dl.assemble(dl.inner(test, trial)*dl.dx)
        test, trial = dl.TestFunction(self.Vh[ADJOINT]), dl.TrialFunction(self.Vh[ADJOINT])
        self.M[ADJOINT] = dl.assemble(dl.inner(test, trial)*dl.dx)
        test, trial = dl.TestFunction(self.Vh[PARAMETER]), dl.TrialFunction(self.Vh[PARAMETER])
        self.M[PARAMETER] = dl.assemble(dl.inner(test, trial)*dl.dx)
        self.Msolver = dl.PETScKrylovSolver('cg', 'jacobi')
        self.Msolver.parameters["maximum_iterations"] = 2000
        self.Msolver.parameters["relative_tolerance"] = 1e-24
        self.Msolver.parameters["absolute_tolerance"] = 1e-24
        self.Msolver.parameters["error_on_nonconvergence"] = True 
        self.Msolver.parameters["nonzero_initial_guess"] = False 
        self.Msolver.set_operator(self.M[PARAMETER])

        self.Prior = jointregularization

        self.alphareg = alphareg

        self.parameters = {'print':False}


    def splitvector(self, x, component="ALL"):
        """ Split parameter x for joint inverse problem into (x1, x2) for each
        independent inverse problem """
        if component == "ALL":
            v1 = [None, None, None]
            v2 = [None, None, None]
            for cc in [STATE, ADJOINT, PARAMETER]:
                v1[cc], v2[cc] = self.splitvector(x[cc], cc)
            return v1, v2
        else:
            fun = vector2Function(x, self.Vh[component])
            fun1, fun2 = fun.split(deepcopy=True)
            return fun1.vector(), fun2.vector()


    def assignvector(self, x1, x2, component="ALL"):
        """ Create a single vector x from (x1, x2) """
        if component == "ALL":
            x = [None, None, None]
            for cc in [STATE, ADJOINT, PARAMETER]:
                x[cc] = self.assignvector(x1[cc], x2[cc], cc)
            return x
        else:
            fun1 = vector2Function(x1, self.model1.problem.Vh[component])
            fun2 = vector2Function(x2, self.model2.problem.Vh[component])
            fun = dl.Function(self.Vh[component])
            dl.assign(fun.sub(0), fun1)
            dl.assign(fun.sub(1), fun2)
            return fun.vector()


    def generate_vector(self, component="ALL"):
        """
        Return the list x=[u,a,p] where:
        - u is any object that describes the state variable
        - a is a Vector object that describes the parameter variable.
          (Need to support linear algebra operations)
        - p is any object that describes the adjoint variable
        
        If component is STATE, PARAMETER, or ADJOINT return x[component]
        """
        if component == "ALL":
            x = [dl.Vector(), dl.Vector(), dl.Vector()]
            self.M[STATE].init_vector(x[STATE], 0)
            self.M[ADJOINT].init_vector(x[ADJOINT], 0)
            self.M[PARAMETER].init_vector(x[PARAMETER], 0)
        else:
            x = dl.Vector()
            self.M[component].init_vector(x, 0)

        return x
            

    def init_parameter(self, a):
        """ Reshape a so that it is compatible with the parameter variable """
        self.M[PARAMETER].init_vector(a, 0)


    def solveFwd(self, out, x, tol=1e-9):
        """ Solve the forward problem  for each inverse problem """
        out1, out2 = self.splitvector(out, STATE)
        x1, x2 = self.splitvector(x, "ALL")

        self.model1.solveFwd(out1, x1, tol)
        self.model2.solveFwd(out2, x2, tol)

        out.zero()
        out.axpy(1.0, self.assignvector(out1, out2, STATE))


    def solveAdj(self, out, x, tol=1e-9):
        """ Solve the adjoint problem for each inverse problem """
        out1, out2 = self.splitvector(out, ADJOINT)
        x1, x2 = self.splitvector(x, "ALL")

        self.model1.solveAdj(out1, x1, tol)
        self.model2.solveAdj(out2, x2, tol)

        out.zero()
        out.axpy(1.0, self.assignvector(out1, out2, ADJOINT))


    def cost(self, x):
        """
        Given the list x = [u,a,p] which describes the state, parameter, and
        adjoint variable compute the cost functional as the sum of 
        the misfit functional and the regularization functional.
        
        Return the list [cost functional, regularization functional, misfit functional]
        
        Note: p is not needed to compute the cost functional
        """        
        x1, x2 = self.splitvector(x, "ALL")

        _, _, misfit1 = self.model1.cost(x1)
        _, _, misfit2 = self.model2.cost(x2)
        misfit = misfit1 + misfit2

        reg = self.Prior.costabvect(x1[PARAMETER], x2[PARAMETER])

        return misfit + self.alphareg*reg, reg, misfit


    def evalGradientParameter(self, x, mg):
        """
        Evaluate the gradient for the variation parameter equation at the point x=[u,a,p].
        Parameters:
        - x = [u,a,p] the point at which to evaluate the gradient.
        - mg the variational gradient (g, atest) being atest a test function in the parameter space
          (Output parameter)
        
        Returns the norm of the gradient in the correct inner product g_norm = sqrt(g,g)
        """ 
        x1, x2 = self.splitvector(x, "ALL")
        mg1, mg2 = self.splitvector(mg, PARAMETER)

        g_n1 = self.model1.evalGradientParameter(x1, mg1)
        g_n2 = self.model2.evalGradientParameter(x2, mg2)
        if self.parameters['print']:
            print '||g1||={}, ||g2||={}'.format(g_n1, g_n2)

        mg.zero()
        mg.axpy(1.0, self.assignvector(mg1, mg2, PARAMETER))
        mg.axpy(self.alphareg, self.Prior.gradabvect(x1[PARAMETER], x2[PARAMETER]))

        g = self.generate_vector(PARAMETER)
        self.Msolver.solve(g, mg)
        g_norm = np.sqrt( g.inner(mg) )

        return g_norm


    def setPointForHessianEvaluations(self, x):
        """
        Specify the point x = [u,a,p] at which the Hessian operator (or the Gauss-Newton approximation)
        need to be evaluated.
        """      
        x1, x2 = self.splitvector(x, "ALL")

        self.model1.setPointForHessianEvaluations(x1)
        self.model2.setPointForHessianEvaluations(x2)

        self.Prior.assemble_hessianab(x1[PARAMETER], x2[PARAMETER])


    def solveFwdIncremental(self, sol, rhs, tol):
        """ Solve the incremental forward problem for a given rhs """
        sol1, sol2 = self.splitvector(sol, STATE)
        rhs1, rhs2 = self.splitvector(rhs, STATE)

        self.model1.solveFwdIncremental(sol1, rhs1, tol)
        self.model2.solveFwdIncremental(sol2, rhs2, tol)

        sol.zero()
        sol.axpy(1.0, self.assignvector(sol1, sol2, STATE))


    def solveAdjIncremental(self, sol, rhs, tol):
        """ Solve the incremental forward problem for a given rhs """
        sol1, sol2 = self.splitvector(sol, ADJOINT)
        rhs1, rhs2 = self.splitvector(rhs, ADJOINT)

        self.model1.solveAdjIncremental(sol1, rhs1, tol)
        self.model2.solveAdjIncremental(sol2, rhs2, tol)

        sol.zero()
        sol.axpy(1.0, self.assignvector(sol1, sol2, ADJOINT))


    def applyC(self, da, out):
        da1, da2 = self.splitvector(da, PARAMETER)
        out1, out2 = self.splitvector(out, STATE)

        self.model1.applyC(da1, out1)
        self.model2.applyC(da2, out2)

        out.zero()
        out.axpy(1.0, self.assignvector(out1, out2, STATE))

    def applyCt(self, dp, out):
        dp1, dp2 = self.splitvector(dp, ADJOINT)
        out1, out2 = self.splitvector(out, PARAMETER)

        self.model1.applyCt(dp1, out1)
        self.model2.applyCt(dp2, out2)

        out.zero()
        out.axpy(1.0, self.assignvector(out1, out2, PARAMETER))

    def applyWuu(self, du, out, gn_approx=False):
        du1, du2 = self.splitvector(du, STATE)
        out1, out2 = self.splitvector(out, ADJOINT)

        self.model1.applyWuu(du1, out1, gn_approx)
        self.model2.applyWuu(du2, out2, gn_approx)

        out.zero()
        out.axpy(1.0, self.assignvector(out1, out2, ADJOINT))

    def applyWua(self, da, out):
        da1, da2 = self.splitvector(da, PARAMETER)
        out1, out2 = self.splitvector(out, ADJOINT)

        self.model1.applyWua(da1, out1)
        self.model2.applyWua(da2, out2)

        out.zero()
        out.axpy(1.0, self.assignvector(out1, out2, ADJOINT))

    def applyWau(self, du, out):
        du1, du2 = self.splitvector(du, STATE)
        out1, out2 = self.splitvector(out, PARAMETER)

        self.model1.applyWau(du1, out1)
        self.model2.applyWau(du2, out2)

        out.zero()
        out.axpy(1.0, self.assignvector(out1, out2, PARAMETER))

    def applyRaa(self, da, out):
        da1, da2 = self.splitvector(da, PARAMETER)
        out1, out2 = self.splitvector(out, PARAMETER)

        self.model1.applyRaa(da1, out1)
        self.model2.applyRaa(da2, out2)

        out.zero()
        out.axpy(1.0, self.assignvector(out1, out2, PARAMETER))


    def applyR(self, da, out):
        da1, da2 = self.splitvector(da, PARAMETER)
        out.zero()
        out.axpy(self.alphareg, self.Prior.hessianab(da1, da2))

        
    def Rsolver(self):        
        return self.Prior.getprecond()


    def mediummisfit(self, a):
        a1, a2 = self.splitvector(a, PARAMETER)
        nd1, nd1p = self.model1.mediummisfit(a1)
        nd2, nd2p = self.model2.mediummisfit(a2)
        return 0.5*(nd1+nd2), 0.5*(nd1p+nd2p)
