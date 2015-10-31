"""
hIPPYlib implements state-of-the-art scalable algorithms for PDE-based
deterministic and Bayesian inverse problems. It builds on FEniCS (a 
parallel finite element element library) [http://fenicsproject.org/]
for the discretization of the PDE and on PETSc [http://www.mcs.anl.gov/petsc/]
for scalable and efficient linear algebra operations and solvers.

For building instructions, see the file INSTALL. Copyright information
and licensing restrictions can be found in the file COPYRIGHT.

The best starting point for new users interested in hIPPYlib's features are the
interactive tutorials in the notebooks folder.

Conceptually, hIPPYlib can be viewed as a toolbox that provides
the building blocks for experimenting new ideas and developing scalable
algorithms for PDE-based deterministic and Bayesian inverse problems.
"""
from variables import *
from modelVerify import *
from NewtonCG import *
from modelTemplate import *
from pd_activeset import *
from assemblePointwiseObservation import assemblePointwiseObservation
from exportPointwiseObservation import exportPointwiseObservation
from linalg import MatMatMult, MatPtAP, to_dense, trace, get_diagonal, estimate_diagonal_inv_coloring, estimate_diagonal_inv2, getColoring, randn_perturb, amg_method
from timeDependentVector import TimeDependentVector
from randomizedEigensolver import singlePass, doublePass, singlePassG, doublePassG
from lowRankOperator import LowRankOperator
from prior import LaplacianPrior, BiLaplacianPrior, ConstrainedBiLaplacianPrior, MollifiedBiLaplacianPrior, LaplaceBeltramiPrior
from posterior import GaussianLRPosterior
from cgsampler import CGSampler
from traceEstimator import TraceEstimator
from expression import code_AnisTensor2D, code_Mollifier
from PDEProblem import PDEProblem, PDEVariationalProblem
from misfit import DistributedStateObservation, PointwiseStateObservation
from model import Model