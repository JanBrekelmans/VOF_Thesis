#pragma once

#include <iostream>

#include "Arrays.h"

struct Grid
{
	Matrix1D x;
	Matrix1D y;
	double dx;
	double dy;
};

void init_grid(Grid& grid, Constants& constants)
{
	double dx = constants.XMAX / constants.NPI;
	double dy = constants.YMAX / constants.NPJ;

	int NPI = constants.NPI;
	int NPJ = constants.NPJ;
	int EX = constants.EX;
	
	grid.dx = dx;
	grid.dy = dy;

	Matrix1D x = double_1D(NPI + 2 * EX);
	Matrix1D y = double_1D(NPJ + 2 * EX);

	for (int I = 0; I < NPI; I++)
	{
		x[I+EX] = 0.5 * dx + I * dx;
	}

	for (int J = 0; J < NPJ; J++)
	{
		y[J+EX] = 0.5 * dy + J * dy;
	}

	for (int m = 0; m < EX; m++)
	{
		x[EX - 1 - m] = x[EX - m] - dx;
		y[EX - 1 - m] = y[EX - m] - dy;

		x[NPI + EX + m] = x[NPI + EX + m - 1] + dx;
		y[NPJ + EX + m] = y[NPJ + EX + m - 1] + dy;
	}

	grid.x = x;
	grid.y = y;
}

void delete_grid(Grid& grid)
{
	delete_Matrix1D(grid.x);
	delete_Matrix1D(grid.y);
}