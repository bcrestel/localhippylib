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

#include "coloring.h"
#include <dolfin/la/PETScVector.h>
#include <dolfin/la/GenericLinearSolver.h>
#include <dolfin/common/Timer.h>

namespace dolfin
{

void cpp_coloring::EstimateDiagonal(const GenericMatrix & A, GenericLinearSolver & Asolver, int k, Vector & diagAinv)
{
	Timer timer("EstimateDiagonal");
	Mat AA = as_type<const PETScMatrix>(A).mat();
	PetscErrorCode ierr;

	Mat AAk;
	Mat AA2;
	timer.start();
	switch(k)
	{
	case 2:
		AAk = AA;
		break;
	case 4:
		::MatMatMult(AA, AA, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &AAk);
		break;
	case 6:
		::MatMatMult(AA, AA, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &AA2);
		::MatMatMult(AA2, AA, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &AAk);
		break;
	case 8:
		::MatMatMult(AA, AA, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &AA2);
		::MatMatMult(AA2, AA2, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &AAk);
		break;
	default:
		std::cout << "Error, only k = 2,4,6,8 supported. Using k = 2\n";
		AAk = AA;
	}

	Vec b, x;
	MatGetVecs(as_type<const PETScMatrix>(A).mat(), &b, &x);

	Vec d_petsc = as_type<PETScVector>(diagAinv).vec();

	ISColoring iscoloring;
    #if PETSC_VERSION_MINOR <= 3
		ierr = MatGetColoring(AAk, MATCOLORINGID, &iscoloring);
    #else
		MatColoring mc;
		MatColoringCreate(AAk,&mc);
		MatColoringSetType(mc,MATCOLORINGID);
		MatColoringSetFromOptions(mc);
		MatColoringApply(mc,&iscoloring);
		MatColoringDestroy(&mc);
    #endif
	timer.stop();
	int size = A.size(0);

	std::cout << "Generate Coloring: " << timer.value() << "\n";
	std::cout << "Number of colors: " << iscoloring->n << "\n";
	std::cout << "Reduction in Solves: " << static_cast<double>(iscoloring->n)/static_cast<double>(size) << "\n";

	timer.start();
	for(int ic = 0; ic < iscoloring->n; ++ic)
	{
		VecSet(b,0.);
		VecSet(x,0.);

		PetscScalar * bb;
		VecGetArray(b,&bb);

		for(int i =0; i < size; ++i)
			if( iscoloring->colors[i] == ic )
				bb[i] = 2.*static_cast<double>(rand() % 2) - 1.;

		VecRestoreArray(b,&bb);
		PETScVector b_fenics(b);
		PETScVector x_fenics(x);
		Asolver.solve(x_fenics, b_fenics);

		PetscScalar * xx;
		PetscScalar * dd;
		VecGetArray(x,&xx);
		VecGetArray(d_petsc,&dd);
		VecGetArray(b,&bb);
		for(int i =0; i < size; ++i)
			if( iscoloring->colors[i] == ic )
			{
				dd[i] = xx[i]*bb[i];
				bb[i] = 0.;
			}

		VecRestoreArray(x,&xx);
		VecRestoreArray(d_petsc,&dd);
		VecRestoreArray(b,&bb);

	}
	timer.stop();
	std::cout << "Solve systems: " << timer.value() << "\n";

/*
	int n = A.size(0);
	std::cout << n << " " << iscoloring->N << "\n";
	std::cout << iscoloring->colors << "\n";
	for(int i =0; i < n; ++i)
		std::cout << iscoloring->colors[i] << "\n";
*/

}


Coloring::Coloring(const GenericMatrix & A, int k)
{

	size = A.size(0);

	Mat AA = as_type<const PETScMatrix>(A).mat();
	PetscErrorCode ierr;

	Mat AAk;
	Mat AA2;
	Mat AA4;

	switch(k)
	{
	case 2:
		AAk = AA;
		break;
	case 4:
		::MatMatMult(AA, AA, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &AAk);
		break;
	case 6:
		::MatMatMult(AA, AA, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &AA2);
		::MatMatMult(AA2, AA, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &AAk);
		MatDestroy(&AA2);
		break;
	case 8:
		::MatMatMult(AA, AA, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &AA2);
		::MatMatMult(AA2, AA2, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &AAk);
		MatDestroy(&AA2);
		break;
	case 12:
		::MatMatMult(AA, AA, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &AA2);
		::MatMatMult(AA2, AA2, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &AA4);
		::MatMatMult(AA2, AA4, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &AAk);
		MatDestroy(&AA2);
		MatDestroy(&AA4);
		break;
	case 16:
		::MatMatMult(AA, AA, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &AA2);
		::MatMatMult(AA2, AA2, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &AA4);
		MatDestroy(&AA2);
		::MatMatMult(AA4, AA4, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &AAk);
		MatDestroy(&AA4);
		break;
	default:
		std::cout << "Error, only k = 2,4,6,8,12,16 supported. Using k = 2\n";
		k = 2;
		AAk = AA;
	}

    #if PETSC_VERSION_MINOR <= 3
		ierr = MatGetColoring(AAk, MATCOLORINGID, &iscoloring);
    #else
		MatColoring mc;
		MatColoringCreate(AAk,&mc);
		MatColoringSetType(mc,MATCOLORINGID);
		MatColoringSetFromOptions(mc);
		MatColoringApply(mc,&iscoloring);
		MatColoringDestroy(&mc);
    #endif
	if(k != 2)
		MatDestroy(&AAk);
}

int Coloring::numberOfColors()
{
	return iscoloring->n;
}
void Coloring::markcolor(int ic, Vector & v, double val)
{
	Vec b = as_type<PETScVector>(v).vec();
	PetscScalar * bb;
	VecGetArray(b,&bb);

	for(int i =0; i < size; ++i)
		if( iscoloring->colors[i] == ic )
			bb[i] = val;

	VecRestoreArray(b,&bb);
}
void Coloring::markcolor_rademaker(int ic, Vector & v)
{
	Vec b = as_type<PETScVector>(v).vec();
	PetscScalar * bb;
	VecGetArray(b,&bb);

	for(int i =0; i < size; ++i)
		if( iscoloring->colors[i] == ic )
			bb[i] = 2.*static_cast<double>(rand() % 2) - 1.;

	VecRestoreArray(b,&bb);
}

void Coloring::copycolor(int ic, const Vector & origin, Vector & dest)
{
	Vec o = as_type<const PETScVector>(origin).vec();
	Vec d = as_type<PETScVector>(dest).vec();
	PetscScalar * oo, *dd;
	VecGetArray(o,&oo);
	VecGetArray(d,&dd);

	for(int i =0; i < size; ++i)
		if( iscoloring->colors[i] == ic )
			dd[i] = oo[i];

	VecRestoreArray(o,&oo);
	VecRestoreArray(d,&dd);
}

void Coloring::x_dot_mult_b(int ic, const Vector & x, Vector & b, Vector & res)
{
	Vec xp = as_type<const PETScVector>(x).vec();
	Vec bp = as_type<PETScVector>(b).vec();
	Vec rp = as_type<PETScVector>(res).vec();
	PetscScalar * xx;
	PetscScalar * bb;
	PetscScalar * rr;
	VecGetArray(xp,&xx);
	VecGetArray(bp,&bb);
	VecGetArray(rp,&rr);
	for(int i =0; i < size; ++i)
		if( iscoloring->colors[i] == ic )
		{
			rr[i] = xx[i]*bb[i];
			bb[i] = 0.;
		}

	VecRestoreArray(xp,&xx);
	VecRestoreArray(bp,&bb);
	VecRestoreArray(rp,&rr);
}

void Coloring::hist(std::vector<int> & h)
{
	h.resize(iscoloring->n);
	std::fill(h.begin(), h.end(), 0);
	for(int i =0; i < size; ++i)
		++h[ iscoloring->colors[i] ];

}

}

