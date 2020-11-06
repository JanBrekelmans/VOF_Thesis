#pragma once

#include <cmath>

#include "Arrays.h"

struct State
{
	Matrix2D u;
	Matrix2D v;
	Matrix2D ut;
	Matrix2D vt;

	Matrix2D p;

	Matrix2D mu;
	Matrix2D rho;

	Matrix2D alpha;
	Matrix2D alphat;

	Matrix2D mxx, mxy, myx, myy;
	Matrix2D nx, ny, k;
};

void init_state(State& state, Grid& grid, Constants& constants)
{
	int NPI = constants.NPI;
	int NPJ = constants.NPJ;
	int EX = constants.EX;

	Matrix1D& x = grid.x;
	Matrix1D& y = grid.y;

	int I, J, i, j;

	state.u = double_2D(NPI + 2 * EX, NPJ + 2 * EX);
	state.v = double_2D(NPI + 2 * EX, NPJ + 2 * EX);
	state.ut = double_2D(NPI + 2 * EX, NPJ + 2 * EX);
	state.vt = double_2D(NPI + 2 * EX, NPJ + 2 * EX);

	state.p = double_2D(NPI + 2 * EX, NPJ + 2 * EX);

	state.mu = double_2D(NPI + 2 * EX, NPJ + 2 * EX);
	state.rho = double_2D(NPI + 2 * EX, NPJ + 2 * EX);

	state.alpha = double_2D(NPI + 2 * EX, NPJ + 2 * EX);
	state.alphat = double_2D(NPI + 2 * EX, NPJ + 2 * EX);

	for (int I = 0; I < NPI + 2 * EX; I++)
	{
		for (int J = 0; J < NPJ + 2 * EX; J++)
		{
			if (y[J] < 0.5 && x[I] < 0.5)
			{
				state.alpha[I][J] = 1.0;
				state.mu[I][J] = constants.MU1;
				state.rho[I][J] = constants.RHO1;
			}
			else
			{
				state.alpha[I][J] = 0.0;
				state.mu[I][J] = constants.MU2;
				state.rho[I][J] = constants.RHO2;
			}
		}
	}

	state.mxx = double_2D(NPI + 2 * EX, NPJ + 2 * EX);
	state.mxy = double_2D(NPI + 2 * EX, NPJ + 2 * EX);
	state.myx = double_2D(NPI + 2 * EX, NPJ + 2 * EX);
	state.myy = double_2D(NPI + 2 * EX, NPJ + 2 * EX);

	state.nx = double_2D(NPI + 2 * EX, NPJ + 2 * EX);
	state.ny= double_2D(NPI + 2 * EX, NPJ + 2 * EX);
	state.k = double_2D(NPI + 2 * EX, NPJ + 2 * EX);
}

void delete_state(State& state)
{
	delete_Matrix2D(state.u);
	delete_Matrix2D(state.v);
	delete_Matrix2D(state.ut);
	delete_Matrix2D(state.vt);

	delete_Matrix2D(state.p);

	delete_Matrix2D(state.mu);
	delete_Matrix2D(state.rho);

	delete_Matrix2D(state.alpha);
	delete_Matrix2D(state.alphat);

	delete_Matrix2D(state.mxx);
	delete_Matrix2D(state.mxy);
	delete_Matrix2D(state.myx);
	delete_Matrix2D(state.myy);

	delete_Matrix2D(state.nx);
	delete_Matrix2D(state.ny);
	delete_Matrix2D(state.k);
}