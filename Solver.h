#pragma once

#include "Arrays.h"

struct Solver_Arrays
{
	Matrix2D aP, aE, aW, aN, aS, b;
};

void init_solver_arrays(Solver_Arrays& solver_arrays, Constants& constants)
{
	int NPI = constants.NPI;
	int NPJ = constants.NPJ;
	int EX = constants.EX;

	solver_arrays.aP = double_2D(NPI + 2 * EX, NPJ + 2 * EX);

	solver_arrays.aE = double_2D(NPI + 2 * EX, NPJ + 2 * EX);
	solver_arrays.aW = double_2D(NPI + 2 * EX, NPJ + 2 * EX);
	solver_arrays.aN = double_2D(NPI + 2 * EX, NPJ + 2 * EX);
	solver_arrays.aS = double_2D(NPI + 2 * EX, NPJ + 2 * EX);

	solver_arrays.b = double_2D(NPI + 2 * EX, NPJ + 2 * EX);
}

void delete_solver_arrays(Solver_Arrays& solver_arrays)
{
	delete_Matrix2D(solver_arrays.aP);

	delete_Matrix2D(solver_arrays.aE);
	delete_Matrix2D(solver_arrays.aW);
	delete_Matrix2D(solver_arrays.aN);
	delete_Matrix2D(solver_arrays.aS);

	delete_Matrix2D(solver_arrays.b);
}