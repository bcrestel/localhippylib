
import dolfin as dl

def targetmediumparameters(Vh, X):
    atrue = dl.interpolate(dl.Expression('\
    0.25*(x[1]>=0.5) + 0.11*(x[1]<0.5)'), Vh)

#    atrue = dl.interpolate(dl.Expression('0.25 - ' + \
#    '(pow(pow(x[0]-0.5,2)+pow(x[1]-0.5,2),0.5)<0.4) * (' + \
#    '0.14*(x[0]<=0.5) + 0.09*(x[0]>0.5) )'), Vh)

#    atrue = dl.interpolate(dl.Expression('2.3 - ' + \
#    '(pow(pow(x[0]-0.5,2)+pow(x[1]-0.5,2),0.5)<0.4) * (' + \
#    '1.6*(x[0]<=0.5) + 0.5*(x[0]>0.5) )'), Vh)
    return atrue,None,None,None,None


def initmediumparameters(Vh, X):
    a0 = dl.interpolate(dl.Constant('0.0'), Vh)
    return a0,None,None,None,None
