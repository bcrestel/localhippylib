"""
Define hippylib model for acoustic wave inverse problem with absorbing boundary
conditions, by wrapping around fenicstools class ObjectiveAcoustic
"""
import dolfin as dl
from variables import STATE, PARAMETER, ADJOINT
from fenicstools.objectiveacoustic import ObjectiveAcoustic


class Acoustic:

    def __init__(self, mpicomm_global, acousticwavePDE, sources, \
    sourcesindex, timestepsindex, regularization=None):
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
        - regularization: object for regularization/prior
    """
    self.objacoustic = ObjectiveAcoustic(mpicomm_global, acousticwavePDE,\
    sources, sourcesindex, timestepsindex, 'a', None)

    self.Prior = regularization

    Vm = self.objacoustic.PDE.Vm
    V = self.objacoustic.PDE.V

    self.problem.Vh[STATE] = V
    self.problem.Vh[ADJOINT] = V
    self.problem.Vh[PARAMETERS] = Vm

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


    def evalGradientParmaeter(self, x, mg, misfit_only=False):
        mga, mgb = self.objacoustic.MG.split(deepcopy=True)
        mg.zero()
        mg.axpy(1.0, mga)

        if not misfit_only:
            mg.axpy(1.0, self.Prior.gradvect(x[PARAMETER]))

        self.Msolver.solve(self.grad, mg)
        mg_norm = np.sqrt(mg.inner(self.grad))
        return mg_norm


    def setPointForHessianEvaluations(self, x):
        a = x[PARAMETER]
        self.objacoustic.update_PDE({'a':a})    #TODO: probably not needed
        self.Prior.assemble_hessian(x[PARAMETER])


    #TODO: def mult, etc....
