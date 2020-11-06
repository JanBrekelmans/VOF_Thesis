#pragma once

#include <cmath>
#include <iostream>
#include <iomanip>

#include "Arrays.h"
#include "Limiter.h"

void bound(State& state, Grid& grid, Constants& constants)
{
	int NPI = constants.NPI;
	int NPJ = constants.NPJ;
	int EX = constants.EX;

	double UNORTH = constants.UNORTH;
	double USOUTH = constants.USOUTH;
	double VEAST = constants.VEAST;
	double VWEST = constants.VWEST;

	Matrix2D& alpha = state.alpha;
	Matrix1D& x = grid.x;
	Matrix1D& y = grid.y;

	for (int I = 0; I < NPI + 2 * EX; I++)
	{
		for (int j = 0; j < EX; j++)
		{
			state.u[I][EX - 1 - j] = 2 * USOUTH - state.u[I][EX + j];
			state.u[I][NPJ + EX + j] = 2 * UNORTH - state.u[I][NPJ + EX - 1 - j];
		}
	}

	for (int J = 0; J < NPJ + 2 * EX; J++)
	{
		for (int i = 0; i < EX; i++)
		{
			state.v[EX - 1 - i][J] = 2 * VEAST - state.v[EX + i][J];
			state.v[NPI + EX + i][J] = 2 * VWEST - state.v[NPI + EX - 1 - i][J];
		}
	}

	for (int I = 0; I < NPI + 2 * EX; I++)
	{
		for (int J = 0; J < EX; J++)
		{
			alpha[I][EX - 1 - J] = alpha[I][EX + J];
			alpha[I][NPJ + EX + J] = alpha[I][NPJ + EX - 1 - J];
		}
	}

	for (int J = 0; J < NPJ + 2 * EX; J++)
	{
		for (int I = 0; I < EX; I++)
		{
			alpha[EX - 1 - I][J] = alpha[EX + I][J];
			alpha[NPI + EX + I][J] = alpha[NPI + EX - 1 - I][J];
		}

	}

}

void calc_ut(State& state, Grid& grid, Constants& constants)
{
	int EX = constants.EX;
	int NPI = constants.NPI;
	int NPJ = constants.NPJ;

	double DT = constants.DT;
	double GX = constants.GX;
	double SIGMA = constants.SIGMA;

	double dx = grid.dx;
	double dy = grid.dy;

	int ISTART = EX + 1;
	int IEND = NPI + EX - 1;
	int JSTART = EX;
	int JEND = NPJ + EX - 1;
	
	int I, J, i, j;
	double A, Ax1, Ax2, Ay1, Ay2;
	double D, Dx1, Dx2, Dy1, Dy2;
	double mu_temp, rho_temp, sigma_temp;

	Matrix2D& u = state.u;
	Matrix2D& v = state.v;
	Matrix2D& ut = state.ut;
	Matrix2D& mu = state.mu;
	Matrix2D& rho = state.rho;
	Matrix2D& alpha = state.alpha;
	Matrix2D& k = state.k;

	for (I = ISTART; I <= IEND; I++)
	{
		i = I;
		for (J = JSTART; J <= JEND; J++)
		{
			j = J;

			Ax1 = (0.25 / dx) * (u[i + 1][J] + u[i][J]) * (u[i + 1][J] + u[i][J]);
			Ax2 = (0.25 / dx) * (u[i][J] + u[i - 1][J]) * (u[i][J] + u[i - 1][J]);

			Ay1 = (0.25 / dy) * (u[i][J + 1] + u[i][J]) * (v[I][j + 1] + v[I - 1][j + 1]);
			Ay2 = (0.25 / dy) * (u[i][J] + u[i][J - 1]) * (v[I][j] + v[I - 1][j]);

			A = Ax1 - Ax2 + Ay1 - Ay2;

			Dx1 = (2 / dx) * mu[I][J] * (u[i + 1][J] - u[i][J]);
			Dx2 = (2 / dx) * mu[I - 1][J] * (u[i][J] - u[i - 1][J]);

			mu_temp = 0.25 * (mu[I][J + 1] + mu[I - 1][J + 1] + mu[I][J] + mu[I - 1][J]);
			Dy1 = mu_temp * ((1 / dy) * (u[i][J + 1] - u[i][J]) + (1 / dx) * (v[I][j + 1] - v[I - 1][j + 1]));
			mu_temp = 0.25 * (mu[I][J] + mu[I - 1][J] + mu[I][J - 1] + mu[I - 1][J - 1]);
			Dy2 = mu_temp * ((1 / dy) * (u[i][J] - u[i][J - 1]) + (1 / dx) * (v[I][j] - v[I - 1][j]));

			D = (1 / dx) * (Dx1 - Dx2) + (1 / dy) * (Dy1 - Dy2);

			rho_temp = 0.5 * (rho[I][J] + rho[I - 1][J]);

			D = D / rho_temp;

			sigma_temp = -0.5 * (1 / dx) * SIGMA * (alpha[I][J] - alpha[I - 1][J]) * (k[I][J] + k[I - 1][J]);

			ut[i][J] = u[i][J] - DT * A + DT * D + DT * GX + DT*sigma_temp;
		}
	}

}

void calc_vt(State& state, Grid& grid, Constants& constants)
{
	int EX = constants.EX;
	int NPI = constants.NPI;
	int NPJ = constants.NPJ;

	double DT = constants.DT;
	double GY = constants.GY;
	double SIGMA = constants.SIGMA;

	double dx = grid.dx;
	double dy = grid.dy;

	int ISTART = EX;
	int IEND = NPI + EX - 1;
	int JSTART = EX + 1;
	int JEND = NPJ + EX - 1;

	int I, J, i, j;
	double A, Ax1, Ax2, Ay1, Ay2;
	double D, Dx1, Dx2, Dy1, Dy2;
	double mu_temp, rho_temp, sigma_temp;

	Matrix2D& u = state.u;
	Matrix2D& v = state.v;
	Matrix2D& vt = state.vt;
	Matrix2D& mu = state.mu;
	Matrix2D& rho = state.rho;
	Matrix2D& alpha = state.alpha;
	Matrix2D& k = state.k;

	for (I = ISTART; I <= IEND; I++)
	{
		i = I;
		for (J = JSTART; J <= JEND; J++)
		{
			j = J;

			Ax1 = (0.25 / dx) * (u[i + 1][J] + u[i + 1][J - 1]) * (v[I + 1][j] + v[I][j]);
			Ax2 = (0.25 / dx) * (u[i][J] + u[i][J - 1]) * (v[I][j] + v[I - 1][j]);

			Ay1 = (0.25 / dy) * (v[I][j + 1] + v[I][j]) * (v[I][j + 1] + v[I][j]);
			Ay2 = (0.25 / dy) * (v[I][j] + v[I][j - 1]) * (v[I][j] + v[I][j - 1]);

			A = Ax1 - Ax2 + Ay1 - Ay2;

			mu_temp = 0.25 * (mu[I + 1][J] + mu[I][J] + mu[I + 1][J - 1] + mu[I][J - 1]);
			Dx1 = mu_temp * ((1 / dy) * (u[i + 1][J] - u[i + 1][J - 1]) + (1 / dx) * (v[I + 1][j] - v[I][j]));
			mu_temp = 0.25 * (mu[I][J] + mu[I - 1][J] + mu[I][J - 1] + mu[I - 1][J - 1]);
			Dx2 = mu_temp * ((1 / dy) * (u[i][J] - u[i][J - 1]) + (1 / dx) * (v[I][j] - v[I - 1][j]));

			Dy1 = (2 / dy) * mu[I][J] * (v[I][j + 1] - v[I][j]);
			Dy2 = (2 / dy) * mu[I][J - 1] * (v[I][j] - v[I][j - 1]);

			D = (1 / dx) * (Dx1 - Dx2) + (1 / dy) * (Dy1 - Dy2);

			rho_temp = 0.5 * (rho[I][J] + rho[I][J - 1]);

			D = D / rho_temp;

			sigma_temp = -0.5 * (1 / dy) * SIGMA * (alpha[I][J] - alpha[I][J - 1]) * (k[I][J] + k[I][J - 1]);

			vt[I][j] = v[I][j] - DT * A + DT * D + DT * GY + DT * sigma_temp;
		}
	}
}

void pcoeffs(State& state, Grid& grid, Solver_Arrays& solver_arrays, Constants& constants)
{
	int EX = constants.EX;
	int NPI = constants.NPI;
	int NPJ = constants.NPJ;

	double DT = constants.DT;

	double dx = grid.dx;
	double dy = grid.dy;

	int ISTART = EX;
	int IEND = NPI + EX - 1;
	int JSTART = EX;
	int JEND = NPJ + EX - 1;

	int I, J, i, j;

	Matrix2D& ut = state.ut;
	Matrix2D& vt = state.vt;
	Matrix2D& rho = state.rho;

	Matrix2D& aP = solver_arrays.aP;
	Matrix2D& aE = solver_arrays.aE;
	Matrix2D& aW = solver_arrays.aW;
	Matrix2D& aN = solver_arrays.aN;
	Matrix2D& aS = solver_arrays.aS;
	Matrix2D& b = solver_arrays.b;

	Matrix2D rt = double_2D(NPI + 2 * EX, NPJ + 2 * EX);

	for (I = 0; I < NPI + 2 * EX; I++)
	{
		for (J = 0; J < NPJ + 2 * EX; J++)
		{
			if ((I < EX) || (I > NPI + EX - 1) || (J < EX) || (J > NPJ + EX - 1))
			{
				rt[I][J] = constants.LARGE;
			}
			else
			{
				rt[I][J] = rho[I][J];
			}
		}
	}


	for (I = ISTART; I <= IEND; I++)
	{
		i = I;
		for (J = JSTART; J <= JEND; J++)
		{
			j = J;
			aP[I][J] = -(1 / (dx * dx)) * (1 / (rt[I + 1][J] + rt[I][J]) + 1 / (rt[I][J] + rt[I - 1][J]))
				- (1 / (dy * dy)) * (1 / (rt[I][J + 1] + rt[I][J]) + 1 / (rt[I][J] + rt[I][J - 1]));

			aE[I][J] = (1 / (dx * dx)) * (1 / (rt[I + 1][J] + rt[I][J]));
			aW[I][J] = (1 / (dx * dx)) * (1 / (rt[I][J] + rt[I - 1][J]));

			aN[I][J] = (1 / (dy * dy)) * (1 / (rt[I][J + 1] + rt[I][J]));
			aS[I][J] = (1 / (dy * dy)) * (1 / (rt[I][J] + rt[I][J - 1]));

			b[I][J] = (0.5 / DT) * ((1 / dx) * (ut[i + 1][J] - ut[i][J]) + (1 / dy) * (vt[I][j + 1] - vt[I][j]));
		}
	}
}

void SOR(Matrix2D& phi, Solver_Arrays& solver_arrays, Constants& constants)
{
	int iter, I, J;
	double max_error = 0.0, current_error;
	double phi_old;

	double OMEGA = constants.OMEGA;
	double MAX_IT = constants.MAX_IT;
	double TOL = constants.TOLERANCE;

	int NPI = constants.NPI;
	int NPJ = constants.NPJ;
	int EX = constants.EX;

	Matrix2D& aP = solver_arrays.aP;

	Matrix2D& aE = solver_arrays.aE;
	Matrix2D& aW = solver_arrays.aW;
	Matrix2D& aN = solver_arrays.aN;
	Matrix2D& aS = solver_arrays.aS;

	Matrix2D& b = solver_arrays.b;

	int ISTART = EX;
	int IEND = NPI + EX - 1;
	int JSTART = EX;
	int JEND = NPJ + EX - 1;

	for (iter = 0; iter < MAX_IT; iter++)
	{
		max_error = 0.0;

		for (I = ISTART; I <= IEND; I++)
		{
			for (J = JSTART; J <= JEND; J++)
			{
				phi_old = phi[I][J];

				phi[I][J] = (1 - OMEGA) * phi[I][J] + (OMEGA / aP[I][J]) * (-aE[I][J] * phi[I + 1][J] - aW[I][J] * phi[I - 1][J]
					- aN[I][J] * phi[I][J + 1] - aS[I][J] * phi[I][J - 1] + b[I][J]);

				current_error = fabs(phi_old - phi[I][J]);

				max_error = current_error > max_error ? current_error : max_error;
			}
		}

		if (max_error < TOL)
		{
			break;
		}
	}

	std::cout << std::setprecision(6);
	std::cout << "Number of iterations : " << std::setw(4) << iter + 1 << "\n";
	std::cout << "Residual : " << std::setw(10) << max_error << "\n";
}

void velocity_correct(State& state, Grid& grid, Constants& constants)
{
	int EX = constants.EX;
	int NPI = constants.NPI;
	int NPJ = constants.NPJ;

	double DT = constants.DT;

	double dx = grid.dx;
	double dy = grid.dy;

	Matrix2D& u = state.u;
	Matrix2D& v = state.v;
	Matrix2D& ut = state.ut;
	Matrix2D& vt = state.vt;
	Matrix2D& mu = state.mu;
	Matrix2D& rho = state.rho;
	Matrix2D& p = state.p;

	int I, J, i, j;
	double rho_temp;

	int ISTART = EX + 1;
	int IEND = NPI + EX - 1;
	int JSTART = EX;
	int JEND = NPJ + EX - 1;

	for (I = ISTART; I <= IEND; I++)
	{
		i = I;
		for (J = JSTART; J <= JEND; J++)
		{
			j = J;

			rho_temp = 0.5 * (rho[I][J] + rho[I - 1][J]);

			u[i][J] = ut[i][J] - (DT / dx) * (p[I][J] - p[I - 1][J]) / rho_temp;
		}
	}

	ISTART = EX;
	IEND = NPI + EX - 1;
	JSTART = EX + 1;
	JEND = NPJ + EX - 1;

	for (I = ISTART; I <= IEND; I++)
	{
		i = I;
		for (J = JSTART; J <= JEND; J++)
		{
			j = J;

			rho_temp = 0.5 * (rho[I][J] + rho[I][J - 1]);

			v[i][J] = vt[i][J] - (DT / dy) * (p[I][J] - p[I][J - 1]) / rho_temp;
		}
	}
}

void calc_normals(State& state, Grid& grid, Constants& constants)
{
	int EX = constants.EX;
	int NPI = constants.NPI;
	int NPJ = constants.NPJ;

	double DT = constants.DT;
	double SMALL = constants.SMALL;

	double dx = grid.dx;
	double dy = grid.dy;

	int ISTART = EX;
	int IEND = NPI + EX - 1;
	int JSTART = EX;
	int JEND = NPJ + EX - 1;

	int I, J, i, j;
	double mtemp;

	Matrix2D& alpha = state.alpha;

	Matrix2D& mxx = state.mxx;
	Matrix2D& mxy = state.mxy;
	Matrix2D& myx = state.myx;
	Matrix2D& myy = state.myy;

	Matrix2D& nx = state.nx;
	Matrix2D& ny = state.ny;

	Matrix2D& k = state.k;

	for (I = ISTART; I <= IEND; I++)
	{
		i = I;
		for (J = JSTART; J <= JEND; J++)
		{
			j = J;
			mxx[i][J] = -(1 / dx) * (alpha[I][J] - alpha[I - 1][J]);
		}
	}

	ISTART = EX;
	IEND = NPI + EX - 1;
	JSTART = EX + 1;
	JEND = NPJ + EX - 1;

	for (I = ISTART; I <= IEND; I++)
	{
		i = I;
		for (J = JSTART; J <= JEND; J++)
		{
			j = J;
			myy[I][j] = -(1 / dy) * (alpha[I][J] - alpha[I][J - 1]);
		}
	}

	ISTART = EX;
	IEND = NPI + EX - 1;
	JSTART = EX;
	JEND = NPJ + EX - 1;

	for (I = ISTART; I <= IEND; I++)
	{
		i = I;
		for (J = JSTART; J <= JEND; J++)
		{
			j = J;
			mxy[i][J] = 0.25 * (myy[I][j + 1] + myy[I - 1][j + 1] + myy[I][j] + myy[I - 1][j]);
		}
	}

	ISTART = EX;
	IEND = NPI + EX - 1;
	JSTART = EX + 1;
	JEND = NPJ + EX - 1;

	for (I = ISTART; I <= IEND; I++)
	{
		i = I;
		for (J = JSTART; J <= JEND; J++)
		{
			j = J;
			myx[I][j] = 0.25 * (mxx[i + 1][J] + mxx[i][J] + mxx[i + 1][J - 1] + mxx[i][J - 1]);
		}
	}

	// Calculate nx
	ISTART = EX;
	IEND = NPI + EX - 1;
	JSTART = EX;
	JEND = NPJ + EX - 1;

	for (I = ISTART; I <= IEND; I++)
	{
		i = I;
		for (J = JSTART; J <= JEND; J++)
		{
			j = J;
			mtemp = mxx[i][J] * mxx[i][J] + mxy[i][J] * mxy[i][J];
			nx[i][J] = mxx[i][J] / (sqrt(mtemp) + SMALL);
		}
	}

	// Calculate ny
	ISTART = EX;
	IEND = NPI + EX - 1;
	JSTART = EX + 1;
	JEND = NPJ + EX - 1;

	for (I = ISTART; I <= IEND; I++)
	{
		i = I;
		for (J = JSTART; J <= JEND; J++)
		{
			j = J;
			
			mtemp = myy[I][j] * myy[I][j] + myx[I][j] * myx[I][j];
			ny[I][j] = myy[I][j] / (sqrt(mtemp) + SMALL);
		}
	}

	// Calculate k
	ISTART = EX;
	IEND = NPI + EX - 1;
	JSTART = EX;
	JEND = NPJ + EX + 1;

	for (I = ISTART; I <= IEND; I++)
	{
		i = I;
		for (J = JSTART; J <= JEND; J++)
		{
			j = J;

			k[I][J] = -(1 / dx) * (nx[i + 1][J] - nx[i][J]) - (1 / dy) * (ny[I][j + 1] - ny[I][j]);
		}
	}
}

void calc_alphat(State& state, Grid& grid, Constants& constants)
{
	int EX = constants.EX;
	int NPI = constants.NPI;
	int NPJ = constants.NPJ;
	
	double DT = constants.DT;
	double SMALL = constants.SMALL;

	double dx = grid.dx;
	double dy = grid.dy;

	Matrix2D& u = state.u;
	Matrix2D& v = state.v;
	Matrix2D& alpha = state.alpha;
	Matrix2D& alphat = state.alphat;

	int ISTART = EX;
	int IEND = NPI + EX - 1;
	int JSTART = EX;
	int JEND = NPJ + EX - 1;

	int I, J, i, j;
	double phiC, phiU, phiD;
	double phi, phi1, phi2;

	for (I = ISTART; I <= IEND; I++)
	{
		i = I;
		for (J = JSTART; J <= JEND; J++)
		{
			j = J;

			alphat[I][J] = 0.0;

			// u velocity
			if (u[i + 1][J] >= 0)
			{
				phiC = alpha[I][J];
				phiU = alpha[I - 1][J];
				phiD = alpha[I + 1][J];
			}
			else
			{
				phiC = alpha[I + 1][J];
				phiU = alpha[I + 2][J];
				phiD = alpha[I][J];
			}
			phi = (phiC - phiU) / (phiD - phiC + SMALL);
			phi1 = phiC + 0.5 * van_Leer(phi) * (phiD - phiC);

			if (u[i][J] >= 0)
			{
				phiC = alpha[I - 1][J];
				phiU = alpha[I - 2][J];
				phiD = alpha[I][J];
			}
			else
			{
				phiC = alpha[I][J];
				phiU = alpha[I + 1][J];
				phiD = alpha[I - 1][J];
			}
			phi = (phiC - phiU) / (phiD - phiC + SMALL);
			phi2 = phiC + 0.5 * van_Leer(phi) * (phiD - phiC);

			alphat[I][J] += (DT / dx) * (phi1 * u[i + 1][J] - phi2 * u[i][J]);

			// v velocity
			if (v[I][j + 1] >= 0)
			{
				phiC = alpha[I][J];
				phiU = alpha[I][J - 1];
				phiD = alpha[I][J + 1];
			}
			else
			{
				phiC = alpha[I][J + 1];
				phiU = alpha[I][J + 2];
				phiD = alpha[I][J];
			}
			phi = (phiC - phiU) / (phiD - phiC + SMALL);
			phi1 = phiC + 0.5 * van_Leer(phi) * (phiD - phiC);

			if (v[I][j] >= 0)
			{
				phiC = alpha[I][J - 1];
				phiU = alpha[I][J - 2];
				phiD = alpha[I][J];
			}
			else
			{
				phiC = alpha[I][J];
				phiU = alpha[I][J + 1];
				phiD = alpha[I][J - 1];
			}
			phi = (phiC - phiU) / (phiD - phiC + SMALL);
			phi2 = phiC + 0.5 * van_Leer(phi) * (phiD - phiC);

			alphat[I][J] += (DT / dy) * (phi1 * v[I][j + 1] - phi2 * v[I][j]);
		}
	}
}

void interface_compression(State& state, Grid& grid, Constants& constants)
{
	int EX = constants.EX;
	int NPI = constants.NPI;
	int NPJ = constants.NPJ;

	double DT = constants.DT;
	double SMALL = constants.SMALL;
	double COMPRESSION = constants.COMPRESSION;

	double dx = grid.dx;
	double dy = grid.dy;

	Matrix2D& u = state.u;
	Matrix2D& v = state.v;
	Matrix2D& alpha = state.alpha;
	Matrix2D& alphat = state.alphat;
	Matrix2D& nx = state.nx;
	Matrix2D& ny = state.ny;

	int ISTART = EX;
	int IEND = NPI + EX - 1;
	int JSTART = EX;
	int JEND = NPJ + EX - 1;

	int I, J, i, j;

	double ax1, ax2, ay1, ay2;
	double alpha_temp;

	for (I = ISTART; I <= IEND; I++)
	{
		i = I;
		for (J = JSTART; J <= JEND; J++)
		{
			j = J;

			if (u[i + 1][J] >= 0)
			{
				alpha_temp = alpha[I][J];
			}
			else
			{
				alpha_temp = alpha[I + 1][J];
			}
			ax1 = -fabs(u[i + 1][J]) * nx[i + 1][J] * alpha_temp * (1 - alpha_temp);

			if (u[i][J] < 0)
			{
				alpha_temp = alpha[I][J];
			}
			else
			{
				alpha_temp = alpha[I - 1][J];
			}
			ax2 = -fabs(u[i][J]) * nx[i][J] * alpha_temp * (1 - alpha_temp);

			if (v[I][j + 1] >= 0)
			{
				alpha_temp = alpha[I][J];
			}
			else
			{
				alpha_temp = alpha[I][J + 1];
			}
			ay1 = -fabs(v[I][j + 1]) * ny[I][j + 1] * alpha_temp * (1 - alpha_temp);

			if (v[I][j] <  0)
			{
				alpha_temp = alpha[I][J];
			}
			else
			{
				alpha_temp = alpha[I][J - 1];
			}
			ay2 = -fabs(v[I][j]) * ny[I][j] * alpha_temp * (1 - alpha_temp);

			alphat[I][J] += COMPRESSION * ((DT / dx) * (ax1 - ax2) + (DT / dy) * (ay1 - ay2));

		}
	}
}

void interface_compression_limited(State& state, Grid& grid, Constants& constants)
{
	int EX = constants.EX;
	int NPI = constants.NPI;
	int NPJ = constants.NPJ;

	double DT = constants.DT;
	double SMALL = constants.SMALL;
	double COMPRESSION = constants.COMPRESSION;

	double dx = grid.dx;
	double dy = grid.dy;

	Matrix2D& u = state.u;
	Matrix2D& v = state.v;
	Matrix2D& alpha = state.alpha;
	Matrix2D& alphat = state.alphat;
	Matrix2D& nx = state.nx;
	Matrix2D& ny = state.ny;

	int ISTART = EX;
	int IEND = NPI + EX - 1;
	int JSTART = EX;
	int JEND = NPJ + EX - 1;

	int I, J, i, j;

	double ax1, ax2, ay1, ay2;
	double alpha_temp, alpha_temp1, alpha_temp2;

	for (I = ISTART; I <= IEND; I++)
	{
		i = I;
		for (J = JSTART; J <= JEND; J++)
		{
			j = J;

			if (nx[i + 1][J] >= 0)
			{
				alpha_temp1 = alpha[I][J];
				alpha_temp2 = alpha[I + 1][J];
			}
			else
			{
				alpha_temp1 = alpha[I + 1][J];
				alpha_temp2 = alpha[I][J];
			}
			alpha_temp = alpha_temp1 + 0.5 * alpha_limiter(alpha_temp1, alpha_temp2) * (alpha_temp2 - alpha_temp1);
			ax1 = -fabs(u[i + 1][J]) * nx[i + 1][J] * alpha_temp * (1 - alpha_temp);

			if (nx[i][J] < 0)
			{
				alpha_temp1 = alpha[I][J];
				alpha_temp2 = alpha[I - 1][J];
			}
			else
			{
				alpha_temp1 = alpha[I - 1][J];
				alpha_temp2 = alpha[I][J];
			}
			alpha_temp = alpha_temp1 + 0.5 * alpha_limiter(alpha_temp1, alpha_temp2) * (alpha_temp2 - alpha_temp1);
			ax2 = -fabs(u[i][J]) * nx[i][J] * alpha_temp * (1 - alpha_temp);

			if (ny[I][j + 1] >= 0)
			{
				alpha_temp1 = alpha[I][J];
				alpha_temp2 = alpha[I][J + 1];
			}
			else
			{
				alpha_temp1 = alpha[I][J + 1];
				alpha_temp2 = alpha[I][J];
			}
			alpha_temp = alpha_temp1 + 0.5 * alpha_limiter(alpha_temp1, alpha_temp2) * (alpha_temp2 - alpha_temp1);
			ay1 = -fabs(v[I][j + 1]) * ny[I][j + 1] * alpha_temp * (1 - alpha_temp);

			if (ny[I][j] < 0)
			{
				alpha_temp1 = alpha[I][J];
				alpha_temp2 = alpha[I][J - 1];
			}
			else
			{
				alpha_temp2 = alpha[I][J - 1];
				alpha_temp1 = alpha[I][J];
			}
			alpha_temp = alpha_temp1 + 0.5 * alpha_limiter(alpha_temp1, alpha_temp2) * (alpha_temp2 - alpha_temp1);
			ay2 = -fabs(v[I][j]) * ny[I][j] * alpha_temp * (1 - alpha_temp);

			alphat[I][J] += COMPRESSION * ((DT / dx) * (ax1 - ax2) + (DT / dy) * (ay1 - ay2));

		}
	}
}

void update_states(State& state, Grid& grid, Constants& constants)
{
	int EX = constants.EX;
	int NPI = constants.NPI;
	int NPJ = constants.NPJ;

	double RHO1 = constants.RHO1;
	double RHO2 = constants.RHO2;

	double MU1 = constants.MU1;
	double MU2 = constants.MU2;

	Matrix2D& alpha = state.alpha;
	Matrix2D& alphat = state.alphat;
	Matrix2D& mu = state.mu;
	Matrix2D& rho = state.rho;

	int I, J;

	for (I = 0; I < NPI + 2 * EX; I++)
	{
		for (J = 0; J < NPJ + 2 * EX; J++)
		{
			alpha[I][J] = alpha[I][J] - alphat[I][J];

			mu[I][J] = MU1 * alpha[I][J] + MU2 * (1 - alpha[I][J]);
			rho[I][J] = RHO1 * alpha[I][J] + RHO2 * (1 - alpha[I][J]);
		}
	}
}