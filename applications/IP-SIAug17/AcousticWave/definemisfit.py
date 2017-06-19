"""
Define misfit function
"""

from hippylib import PointwiseStateObservation
import numpy as np

def defmisfit(Vh, STATE):
    nbobsperdir = 50
    targets = np.array([ [float(i)/(nbobsperdir+1), float(j)/(nbobsperdir+1)] \
    for i in range(1, nbobsperdir+1) for j in range(1, nbobsperdir+1)])
    misfit = PointwiseStateObservation(Vh[STATE], targets)
    return misfit, targets
