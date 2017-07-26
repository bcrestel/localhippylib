import sys
import numpy as np
import dolfin as dl

from hippylib import vector2Function
from hippylib import STATE, ADJOINT, PARAMETER
from hippylib import ReducedHessian

from fenicstools.linalg.splitandassign import SplitAndAssign, SplitAndAssigni
from fenicstools.miscfenics import createMixedFS, createMixedFSi


class JointModel:
    """ This class builds a model for joint inversion in the case of two
    independent physics connected only via the regularization term (or prior) """

    def __init__(self, model1, model2, jointregularization, \
        alphareg=1.0, parameters=[]):
        """ Provides models model1 and model2 with ZeroPrior as their
        regularization, and provide jointregularization. """

        self.model1 = model1
        self.model2 = model2

        # individual models shall not have their own regularization!
        assert self.model1.Prior.isZeroPrior()
        assert self.model2.Prior.isZeroPrior()

        self.Vh = [None, None, None]
        self.Vh[STATE] = createMixedFS(self.model1.problem.Vh[STATE], self.model2.problem.Vh[STATE])
        self.Vh[ADJOINT] = createMixedFS(self.model1.problem.Vh[ADJOINT], self.model2.problem.Vh[ADJOINT])
        self.Vh[PARAMETER] = createMixedFS(self.model1.problem.Vh[PARAMETER], self.model2.problem.Vh[PARAMETER])

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

        self.GN = False
        self.tol = 1e-9

        self.Prior = jointregularization
        self.alphareg = alphareg

        self.parameters = {'print':False, 'splitassign':True}
        self.parameters.update(parameters)

        if self.parameters['splitassign']:
            try:
                self.splitassign = [None, None, None]
                for cc in [STATE, ADJOINT, PARAMETER]:
                    self.splitassign[cc] = SplitAndAssign( \
                    self.model1.problem.Vh[cc], self.model2.problem.Vh[cc], \
                    self.M[cc].mpi_comm())
                if self.parameters['print']:
                    print '[JointModel] Using SplitAndAssign'
            except:
                self.parameters['splitassign'] = False
                if self.parameters['print']:
                    print '[JointModel] NOT using SplitAndAssign'
        else:
            if self.parameters['print']:
                print '[JointModel] NOT using SplitAndAssign'


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
            if self.parameters['splitassign'] == True:
                return self.splitassign[component].split(x)
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
            if self.parameters['splitassign'] == True:
                return self.splitassign[component].assign(x1, x2)
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
        if self.parameters['splitassign']:
            mg1 = self.model1.generate_vector(PARAMETER)
            mg2 = self.model2.generate_vector(PARAMETER)
        else:
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


    def mult(self, x, y):
        """
        Compute Hessian-vector product with x, and save it to y
        """
        x1, x2 = self.splitvector(x, PARAMETER)
        y1, y2 = self.splitvector(y, PARAMETER)

        # model 1
        H1 = ReducedHessian(self.model1, self.tol, self.GN, True)
        H1.mult(x1, y1)
        # model 2
        H2 = ReducedHessian(self.model2, self.tol, self.GN, True)
        H2.mult(x2, y2)
        # combine them
        y.zero()
        y.axpy(1.0, self.assignvector(y1, y2, PARAMETER))


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

    def getPDEcounts(self):
        return self.model1.getPDEcounts() + self.model2.getPDEcounts()





class JointModeli:
    """ This class builds a model for joint inversion in the case of n
    independent physics connected only via the regularization term (or prior) """

    #TODO: modify jointregularization to accept single vector (no matter nb of param)
    def __init__(self, models, jointregularization, \
        alphareg=1.0, parameters=[]):
        """ Provides list of models with ZeroPrior as their
        regularization, and provide jointregularization. """

        self.models = models

        # individual models shall not have their own regularization!
        for model in self.models:
            assert self.model.Prior.isZeroPrior()

        Vhstate, Vhadj, Vhparam = [], [], []
        for model in self.models:
            Vhstate.append(model.problem.Vh[STATE])
            Vhadj.append(model.problem.Vh[ADJOINT])
            Vhparam.append(model.problem.Vh[PARAMETER])
        Vhs = [None, None, None]
        Vhs[STATE] = Vhstate
        Vhs[ADJOINT] = Vhadj
        Vhs[PARAMETER] = Vhparam

        self.Vh = [None, None, None]
        for cc in [STATE, ADJOINT, PARAMETER]:
            self.Vh[cc] = createMixedFSi(Vhs[cc])

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

        self.GN = False
        self.tol = 1e-9

        self.Prior = jointregularization
        self.alphareg = alphareg

        self.parameters = {'print':False, 'splitassign':True}
        self.parameters.update(parameters)

        if self.parameters['splitassign']:
            try:
                self.splitassign = [None, None, None]
                for cc in [STATE, ADJOINT, PARAMETER]:
                    self.splitassign[cc] = SplitAndAssigni(Vhs[cc])
                if self.parameters['print']:
                    print '[JointModel] Using SplitAndAssigni'
            except:
                self.parameters['splitassign'] = False
                if self.parameters['print']:
                    print '[JointModel] NOT using SplitAndAssign'
        else:
            if self.parameters['print']:
                print '[JointModel] NOT using SplitAndAssign'


    def splitvector(self, x, component="ALL"):
        """ Split parameter x for joint inverse problem into (x1, ..., xn) for each
        independent inverse problem """
        if component == "ALL":
            vv = []
            for mm in self.models:
                vv.append([None, None, None])
            for cc in [STATE, ADJOINT, PARAMETER]:
                vcc = self.splitvector(x[cc], cc)
                for ii, vvcc in enumerate(vcc):
                    vv[ii][cc] = vvcc
            return vv
        else:
            if self.parameters['splitassign'] == True:
                return self.splitassign[component].split(x)
            else:
                sys.exit(1)
#                fun = vector2Function(x, self.Vh[component])
#                funs = fun.split(deepcopy=True)
#                funsout = []
#                for ff in funs:
#                    funsout.append(funs.vector())
#                return funsout


    def assignvector(self, xs, component="ALL"):
        """ Create a single vector x from (x1,..,xn) """
        if component == "ALL":
            x = [None, None, None]
            for cc in [STATE, ADJOINT, PARAMETER]:
                xscc = []
                for xsii in xs:
                    xscc.append(xsii[cc])
                x[cc] = self.assignvector(xscc, cc)
            return x
        else:
            if self.parameters['splitassign'] == True:
                return self.splitassign[component].assign(xs)
            else:
                sys.exit(1)
#                fun1 = vector2Function(x1, self.model1.problem.Vh[component])
#                fun2 = vector2Function(x2, self.model2.problem.Vh[component])
#                fun = dl.Function(self.Vh[component])
#                dl.assign(fun.sub(0), fun1)
#                dl.assign(fun.sub(1), fun2)
#                return fun.vector()


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
        outs = self.splitvector(out, STATE)
        xs = self.splitvector(x, "ALL")

        for modelii, outii, xsii in zip(self.models, outs, xs):
            self.modelii.solveFwd(outii, xsii, tol)
    
        out.zero()
        out.axpy(1.0, self.assignvector(outs, STATE))


    def solveAdj(self, out, x, tol=1e-9):
        """ Solve the adjoint problem for each inverse problem """
        outs = self.splitvector(out, ADJOINT)
        xs = self.splitvector(x, "ALL")

        for modelii, outii, xsii in zip(self.models, outs, xs):
            self.modelii.solveAdj(outii, xsii, tol)
    
        out.zero()
        out.axpy(1.0, self.assignvector(outs, ADJOINT))


    def cost(self, x):
        """
        Given the list x = [u,a,p] which describes the state, parameter, and
        adjoint variable compute the cost functional as the sum of 
        the misfit functional and the regularization functional.
        
        Return the list [cost functional, regularization functional, misfit functional]
        
        Note: p is not needed to compute the cost functional
        """        
        xs = self.splitvector(x, "ALL")
        misfit = 0.0

        for modelii, xsii in zip(self.models, xs):
            _, _, misfitii = self.modelii.cost(xsii)
            misfit += misfitii

        reg = self.Prior.costabvect(x[PARAMETER])

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
        xs = self.splitvector(x, "ALL")
        mgs = []
        for modelii, xsii in zip(self.models, xs):
            mgii = modelii.generate_vector(PARAMETER)
            g_nii = modelii.evalGradientParameter(xsii, mgii)
            if self.parameters['print']:
                print '||g||={},'.format(g_nii),
            mgs.append(mgii)
        if self.parameters['print']:    print ''

        mg.zero()
        mg.axpy(1.0, self.assignvector(mgs, PARAMETER))
        mg.axpy(self.alphareg, self.Prior.gradabvect(xs[PARAMETER])

        g = self.generate_vector(PARAMETER)
        self.Msolver.solve(g, mg)
        g_norm = np.sqrt( g.inner(mg) )

        return g_norm


    def setPointForHessianEvaluations(self, x):
        """
        Specify the point x = [u,a,p] at which the Hessian operator (or the Gauss-Newton approximation)
        need to be evaluated.
        """      
        xs = self.splitvector(x, "ALL")
        for modelii, xsii in zip(self.models, xs):
            modelii.setPointForHessianEvaluations(xsii)

        self.Prior.assemble_hessianab(xs[PARAMETER])


    def mult(self, x, y):
        """
        Compute Hessian-vector product with x, and save it to y
        """
        xs = self.splitvector(x, PARAMETER)
        ys = self.splitvector(y, PARAMETER)

        for modelii, xsii, ysii in zip(self.models, xs, ys):
            Hii = ReducedHessian(modelii, self.tol, self.GN, True)
            #TODO: check that ysii is updated
            Hii.mult(xsii, ysii)

        y.zero()
        y.axpy(1.0, self.assignvector(ys, PARAMETER))


    def applyR(self, da, out):
        out.zero()
        out.axpy(self.alphareg, self.Prior.hessianab(da))

        
    def Rsolver(self):        
        return self.Prior.getprecond()


    def mediummisfit(self, a):
        a_s = self.splitvector(a, PARAMETER)
        nd = 0.0
        ndp = 0.0
        for modelii, a_sii in zip(self.models, a_s):
            ndii, ndiip = modelii.mediummisfit(a_sii)
            nd += ndii
            ndp += ndiip
        return 0.5*nd, 0.5*ndp


    def getPDEcounts(self):
        PDEc = 0.0
        for modelii in self.models:
            PDEc += modelii.getPDEcounts()
        return PDEc
