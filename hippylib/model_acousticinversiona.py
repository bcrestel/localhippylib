"""
Define hippylib model for acoustic wave inverse problem with absorbing boundary
conditions, by wrapping around fenicstools class ObjectiveAcoustic
"""
import numpy as np
import dolfin as dl

from variables import STATE, PARAMETER, ADJOINT

from fenicstools.objectiveacoustic import ObjectiveAcoustic
from fenicstools.miscfenics import createMixedFS



class FakeProblem:
    """
    Implement fake 'problem' class that only has Vh as a member
    """
    def __init__(self):
        self.Vh = [None, None, None]


class ModelAcoustic:

    def __init__(self, mpicomm_global, acousticwavePDE, sources, sourcesindex,\
    timestepsindex, obsop, regularization, atrue=None, btrue=None):
        """
        Arguments:
            - mpicomm_global = MPI communicator to average gradient and Hessian-vect
            contributions
            - acousticwavePDE = object from fenicstools' class
            acousticwave.AcousticWave
            - sources = seismic source for PDE
            - sourcesindex = list of source numbers to be run by this MPI proc
            - timestepsindex = range of time steps to be computed by this MPI proc
              (for gradient and Hessian-vect)
            - obsop = observation operator
            - atrue: target medium for parameter a (Function)
            - btrue: target medium for parameter b (Function)
            - regularization: object for regularization/prior
        """
        self.objacoustic = ObjectiveAcoustic(mpicomm_global, acousticwavePDE,\
        sources, sourcesindex, timestepsindex, 'a', None)
        self.objacoustic.alpha_reg = 0.0    # belt AND hangers
        self.objacoustic.obsop = obsop
        self.GN = False
        self.tol = 1e-9

        self.Prior = regularization

        Vm = self.objacoustic.PDE.Vm
        V = self.objacoustic.PDE.V
        VmVm = createMixedFS(Vm, Vm)

        self.problem = FakeProblem()
        self.problem.Vh[STATE] = V
        self.problem.Vh[ADJOINT] = V
        self.problem.Vh[PARAMETER] = Vm

        self.M = [None, None, None]
        test, trial = dl.TestFunction(V), dl.TrialFunction(V)
        self.M[STATE] = dl.assemble(dl.inner(test, trial)*dl.dx)
        test, trial = dl.TestFunction(V), dl.TrialFunction(V)
        self.M[ADJOINT] = dl.assemble(dl.inner(test, trial)*dl.dx)
        test, trial = dl.TestFunction(Vm), dl.TrialFunction(Vm)
        self.M[PARAMETER] = dl.assemble(dl.inner(test, trial)*dl.dx)
        self.Msolver = dl.PETScKrylovSolver('cg', 'jacobi')
        self.Msolver.parameters["maximum_iterations"] = 2000
        self.Msolver.parameters["relative_tolerance"] = 1e-24
        self.Msolver.parameters["absolute_tolerance"] = 1e-24
        self.Msolver.parameters["error_on_nonconvergence"] = True 
        self.Msolver.parameters["nonzero_initial_guess"] = False 
        self.Msolver.set_operator(self.M[PARAMETER])

        self.grad = self.generate_vector(PARAMETER)
        self.a = dl.Function(Vm)
        self.x_ab = dl.Function(VmVm)
        self.y_ab = dl.Function(VmVm)

        self.atrue = atrue
        self.btrue = btrue
        self.abtrue = dl.Function(VmVm)
        self.atruen =\
        np.sqrt(self.atrue.vector().inner(self.M[PARAMETER]*self.atrue.vector()))


    def generate_vector(self, component = "ALL"):
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
        self.M[PARAMETER].init_vector(a, 0)


    def cost(self, x):
        cost_misfit = self.objacoustic.cost_misfit
        cost_reg = self.Prior.costvect(x[PARAMETER])
        cost = cost_misfit + cost_reg
        return cost, cost_reg, cost_misfit


    def solveFwd(self, out, x, tol=1e-9):
        out.zero()

        self.objacoustic.update_PDE({'a':x[PARAMETER]})
        self.objacoustic.solvefwd_cost()


    def solveAdj(self, out, x, tol=1e-9):
        out.zero()

        self.objacoustic.solveadj_constructgrad()


    def evalGradientParameter(self, x, mg, misfit_only=False):
        mga, mgb = self.objacoustic.MG.split(deepcopy=True)
        mg.zero()
        mg.axpy(1.0, mga.vector())

        if not misfit_only:
            mg.axpy(1.0, self.Prior.gradvect(x[PARAMETER]))

        self.Msolver.solve(self.grad, mg)
        mg_norm = np.sqrt(mg.inner(self.grad))
        return mg_norm


    def setPointForHessianEvaluations(self, x):
        a = x[PARAMETER]
        self.objacoustic.update_PDE({'a':a})    #TODO: probably not needed
        self.Prior.assemble_hessian(x[PARAMETER])


    def mult(self, x, y):
        """
        Compute Hessian-vect product with x, and return to y
        """
        self.a.vector().zero()
        self.a.vector().axpy(1.0, x)

        self.x_ab.vector().zero()
        dl.assign(self.x_ab.sub(0), self.a)

        self.objacoustic.GN = self.GN
        self.objacoustic.mult(self.x_ab.vector(), self.y_ab.vector())

        ya, yb = self.y_ab.split(deepcopy=True)
        y.zero()
        y.axpy(1.0, ya.vector())


    def applyR(self, da, out):
        """
        Apply the regularization R to a (incremental) parameter variable.
        out = R da
        Parameters:
        - da the (incremental) parameter variable
        - out the action of R on da
        
        Note: this routine assumes that out has the correct shape.
        """
        out.zero()
        out.axpy(1.0, self.Prior.hessian(da))


    def Rsolver(self):
        """
        Return an object Rsovler that is a suitable solver for the regularization
        operator R.
        """
        return self.Prior.getprecond()


    def mediummisfit(self, m):
        if self.atrue == None:
            return -99, -99
        else:
            self.abtrue.vector().zero()
            dl.assign(self.abtrue.sub(0), self.atrue)
            dl.assign(self.abtrue.sub(1), self.btrue)
            mma, mmb = self.objacoustic.mediummisfit(self.abtrue)
            assert mmb < 1e-24, mmb
            return mma, 100.0*mma/self.atruen


    def generate_synthetic_obs(self, SNRdB=0.0):
        """
        Generate synthetic observations at the target parameter
        """
        self.objacoustic.update_PDE({'a':self.atrue, 'b':self.btrue})
        self.objacoustic.solvefwd()
        DD = self.objacoustic.Bp[:]
        if SNRdB > 0.0:
            np.random.seed(11)
            # Generate random components for all src (even if not owned)
            RAND = []
            nbobspt = self.objacoustic.obsop.PtwiseObs.nbPts
            nbt = self.objacoustic.PDE.Nt + 1
            Pt = self.objacoustic.fwdsource[1]
            for ii in range(len(Pt.src_loc)):
                RAND.append(np.random.randn(nbobspt*nbt).reshape((nbobspt, nbt)))
            sources = self.objacoustic.srcindex
            RAND = RAND[sources[0]:sources[-1]+1]
            mpiglobalrank = dl.MPI.rank(self.objacoustic.mpicomm_global)
            # Add noise
            for ii in range(len(DD)):
                dd = DD[ii]
                rndnoise = RAND[ii]
                mu = np.abs(dd).mean(axis=1)
                sigmas = mu/(10**(SNRdB/10.))
                DD[ii] = dd + sigmas.reshape((nbobspt,1))*rndnoise
                print '[ModelAcoustic] source={}, mpiglobalrank={}, sigmas={}, |rndnoise|={}'.format(\
                ii, mpiglobalrank, sigmas.sum()/len(sigmas), (rndnoise**2).sum().sum())
                dl.MPI.barrier(dl.mpi_comm_world())
        self.objacoustic.dd = DD


    #TODO: create this
    def getPDEcounts(self):
        return 0
