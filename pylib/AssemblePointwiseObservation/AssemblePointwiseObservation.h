#include <dolfin/function/FunctionSpace.h>
#include <dolfin/la/PETScMatrix.h>
#include <dolfin/la/Matrix.h>

namespace dolfin
{

class PointwiseObservation
{
public:
	PointwiseObservation(const FunctionSpace & Vh, const Array<double> & targets);
	Matrix GetMatrix();

private:
	Mat mat;
};

}
