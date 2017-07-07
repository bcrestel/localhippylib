"""
Medium parameters
"""

import dolfin as dl


def createparam(Vl, DD, Vt, Vb, Vp):
    c = dl.interpolate(dl.Expression('\
    (x[0]>=LL)*(x[0]<=RRR)*(x[1]>=BB)*(x[1]<=TT)*vp\
    + (1.0-(x[0]>=LL)*(x[0]<=RRR)*(x[1]>=BB)*(x[1]<=TT))*\
    (vt*(x[1]>=0.5) + vb*(x[1]<0.5))',\
    LL=0.5-0.5*DD, RRR=0.5+0.5*DD, \
    BB=0.5, TT=0.7,\
    vp=Vp, vt=Vt, vb=Vb, degree=10), Vl)
    return c

DD = 0.5
# medium parameters:
CC = [2.0, 3.0, 2.5]
AAa = []
for cc in CC:
    AAa.append(1./cc**2)
fact = 2.0
AAp = [fact*AAa[0], AAa[1], fact*AAa[0]]

def targetmediumparameters(Vl, X=1.0, myplot=None):
    """
    Arguments:
        Vl = function space
        X = 0.0->elliptic, 1.0->acoustic
    """
    if X > 0.5: AA = AAa
    else:   AA = AAp

    af = createparam(Vl, DD, AA[0], AA[1], AA[2])
    if not myplot == None:
        myplot.set_varname('alpha_target_X'+ str(X))
        myplot.plot_vtk(af)
    
    bf = dl.interpolate(dl.Constant('1.0'), Vl)

    return af, bf


def initmediumparameters(Vl, X=1.0, myplot=None):
    if X > 0.5: AA = AAa
    else:   AA = AAp

    a0 = dl.interpolate(dl.Constant(str(AA[0])), Vl)
    if not myplot == None:
        myplot.set_varname('alpha_init_X'+ str(X))
        myplot.plot_vtk(a0)
    b0 = dl.interpolate(dl.Constant('1.0'), Vl)

    return a0, b0



def loadparameters(LARGE):
    if LARGE:
        Nxy = 20
        Dt = 1.0e-3
        fpeak = 4.0
        t0, t1, t2, tf = 0.0, 0.1, 0.9, 1.0
    else:
        Nxy = 20
        Dt = 1.5e-3
        fpeak = 2.0
        t0, t1, t2, tf = 0.0, 0.1, 1.7, 1.8
    return Nxy, Dt, fpeak, t0, t1, t2, tf
