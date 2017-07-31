"""
Define misfit function
"""

from hippylib import PointwiseStateObservation
import numpy as np

def defmisfit(Vh, STATE):
    nbobsperdir=50
    targets1 = np.array([ [float(i)/(nbobsperdir+1), float(j)/(nbobsperdir+1)] \
    for i in range((nbobsperdir+2)/2, nbobsperdir+1) \
    for j in range((nbobsperdir+2)/2, nbobsperdir+1)])
    targets2 = np.array([ [float(i)/(nbobsperdir+1), float(j)/(nbobsperdir+1)] \
    for i in range(1, nbobsperdir+1) for j in range(1, nbobsperdir+1)])
    misfit1 = PointwiseStateObservation(Vh[STATE], targets1)
    misfit2 = PointwiseStateObservation(Vh[STATE], targets2)

    return misfit1, misfit2, targets1, targets2
