import dolfin as dl
import sys
sys.path.append( "../" )
from pylib import *
import numpy as np

code_AnisTensor2D = '''
class AnisTensor2D : public Expression
{
public:

  AnisTensor2D() :
  Expression(2,2),
  theta0(1.),
  theta1(1.),
  alpha(0)
  {

  }

void eval(Array<double>& values, const Array<double>& x) const
  {
     double sa = sin(alpha);
     double ca = cos(alpha);
     double c00 = theta0*sa*sa + theta1*ca*ca;
     double c01 = (theta0 - theta1)*sa*ca;
     double c11 = theta0*ca*ca + theta1*sa*sa;
  
     values[0] = c00;
     values[1] = c01;
     values[2] = c01;
     values[3] = c11;
  }
  
  double theta0;
  double theta1;
  double alpha;
  
};
'''

code_Mollifier = '''
class Mollifier : public Expression
{

public:

  Mollifier() :
  Expression(),
  nlocations(0),
  locations(nlocations),
  l(1),
  o(2),
  theta0(1),
  theta1(1),
  alpha(0)
  {
  }

void eval(Array<double>& values, const Array<double>& x) const
  {
        double sa = sin(alpha);
        double ca = cos(alpha);
        double c00 = theta0*sa*sa + theta1*ca*ca;
        double c01 = (theta0 - theta1)*sa*ca;
        double c11 = theta0*ca*ca + theta1*sa*sa;
        
        int ndim(2);
        Array<double> dx(ndim);
        double e(0), val(0);
        for(int ip = 0; ip < nlocations; ++ip)
        {
            for(int idim = 0; idim < ndim; ++idim)
                dx[idim] = x[idim] - locations[2*ip+idim];
                
            e = dx[0]*dx[0]*c00 + dx[1]*dx[1]*c11 + 2*dx[0]*dx[1]*c01;
            val += exp( -pow(e/(l*l), .5*o) );
        }
        values[0] = val;
  }
  
  void addLocation(double x, double y) { locations.push_back(x); locations.push_back(y); ++nlocations;}
  
  double l;
  double o;
  
  double theta0;
  double theta1;
  double alpha;
  
  private:
    int nlocations;
    std::vector<double> locations;
  
};
'''

if __name__ == "__main__":
    dl.set_log_active(False)
    ndim = 2
    nx = 64
    ny = 64
    mesh = dl.UnitSquareMesh(nx, ny)
    Vh = dl.FunctionSpace(mesh, "CG", 1)
    
    T = dl.Expression(code_AnisTensor2D)
    T.theta0 = 2.
    T.theta1 = .5
    T.alpha = math.pi/4
    
    u = dl.TrialFunction(Vh)
    v = dl.TestFunction(Vh)
    a = dl.inner(T*dl.grad(u), dl.grad(v))*dl.dx
    
    A = dl.assemble(a)
    
    e = dl.Expression(code_Mollifier)
    e.o = 2
    e.l = 0.2
    e.addLocation(.1,.1)
    e.addLocation(.1,.9)
    e.addLocation(.5,.5)
    e.addLocation(.9,.1)
    e.addLocation(.9,.9)
    e.theta0 = 1./T.theta0
    e.theta1 = 1./T.theta1
    e.alpha  = T.alpha
    
    m = dl.interpolate(e, Vh)
    
    dl.plot(m)
    dl.interactive()