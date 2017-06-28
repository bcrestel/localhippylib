"""
Target medium for test case "coincide1"
"""

import dolfin as dl

def targetmedium(Vh, PARAMETER):
    a1true = dl.interpolate(dl.Expression('log(10 - ' + \
    '(pow(pow(x[0]-0.5,2)+pow(x[1]-0.5,2),0.5)<0.4) * 8 )', degree=10), Vh[PARAMETER])
    a2true = dl.interpolate(dl.Expression('log(10 - ' + \
    '(pow(pow(x[0]-0.5,2)+pow(x[1]-0.5,2),0.5)<0.4) * (' + \
    '8*(x[0]<=0.5) + 4*(x[0]>0.5) ))', degree=10), Vh[PARAMETER])

    return a1true, a2true
