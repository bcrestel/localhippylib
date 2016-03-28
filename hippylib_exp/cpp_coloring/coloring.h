/* Copyright (c) 2016, The University of Texas at Austin & University of
 * California, Merced.
 *
 * All Rights reserved.
 * See file COPYRIGHT for details.
 *
 * This file is part of the hIPPYlib library. For more information and source
 * code availability see https://hippylib.github.io.
 *
 * hIPPYlib is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License (as published by the Free
 * Software Foundation) version 3.0 dated June 2007.
*/

#include <dolfin/la/GenericMatrix.h>
#include <dolfin/la/Matrix.h>
#include <dolfin/la/Vector.h>
#include <dolfin/la/PETScMatrix.h>

namespace dolfin
{

	class cpp_coloring
	{
	public:
		void EstimateDiagonal(const GenericMatrix & A, GenericLinearSolver & Asolver, int k, Vector & diagAinv);
	};

	class Coloring
	{
	public:
		Coloring(const GenericMatrix & A, int k);
		int numberOfColors();
		void markcolor(int color, Vector & v, double val);
		void markcolor_rademaker(int color, Vector & v);
		void copycolor(int color, const Vector & origin, Vector & dest);
		void x_dot_mult_b(int color, const Vector & x, Vector & b, Vector & res);
		void hist(std::vector<int> & h);
	private:
		ISColoring iscoloring;
		int size;
	};

}
