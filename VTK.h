#pragma once

#include <iostream>
#include <fstream>
#include <sstream>

void write_to_VTK(State& state, Grid& grid, Constants& constants, int iter, double time)
{
	int EX = constants.EX;
	int NPI = constants.NPI;
	int NPJ = constants.NPJ;

	double SMALLEST_WRITE_VALUE = constants.SMALLEST_WRITE_VALUE;
	double LARGEST_WRITE_VALUE = constants.LARGEST_WRITE_VALUE;

	double dx = grid.dx;
	double dy = grid.dy;

	int I, J, i, j;

	double value;
	int sign;

	Matrix2D& u = state.u;
	Matrix2D& v = state.v;
	Matrix2D& p = state.p;
	Matrix2D& alpha = state.alpha;
	Matrix2D& mu = state.mu;
	Matrix2D& rho = state.rho;

	std::ofstream file;

	std::string file_name = "Iteration" + std::to_string(iter) + ".vtk";
	std::string path_name = "VTK/" + file_name;

	file.open(path_name, std::ios::out);

	file << std::setprecision(5);

	file << "# vtk DataFile Version 2.0" << std::endl
		<< "Results" << std::endl
		<< "ASCII" << std::endl;
	file << "DATASET STRUCTURED_POINTS" << std::endl
		<< "DIMENSIONS " << NPI + 1 << " " << NPJ + 1 << " " << 1 << std::endl
		<< "ORIGIN 0.0 0.0 0.0" << std::endl
		<< "SPACING " << std::to_string(dx) << " " << std::to_string(dy) << " " << std::to_string(1.0) << std::endl << std::endl;

	file << "FIELD FieldData 1 \n";
	file << "TIME 1 1 double\n"
		<< time << "\n";

	file << "CELL_DATA " << NPI * NPJ << std::endl;
	file << "SCALARS pressure float 1" << std::endl;
	file << "LOOKUP_TABLE default" << std::endl;
	
	int ISTART = EX;
	int IEND = NPI + EX - 1;
	int JSTART = EX;
	int JEND = NPJ + EX - 1;

	for (J = JSTART; J <= JEND; J++)
	{
		j = J;;
		for (I = ISTART; I <= IEND; I++)
		{
			i = I;

			sign = p[I][J] >= 0 ? 1 : -1;
			value = fabs(p[I][J]) < SMALLEST_WRITE_VALUE ? 0 : p[I][J];
			value = fabs(p[I][J]) > LARGEST_WRITE_VALUE ? sign * LARGEST_WRITE_VALUE : p[I][J];

			file << value << "\t";
		}
		file << "\n";
	}

	file << "\n";

	file << "VECTORS U float\n";
	for (J = JSTART; J <= JEND; J++)
	{
		j = J;
		for (I = ISTART; I <= IEND; I++)
		{
			i = I;

			value = 0.5 * (u[i + 1][J] + u[i][J]);
			sign = value >= 0 ? 1 : -1;
			value = fabs(value) < SMALLEST_WRITE_VALUE ? 0 : value;
			value = fabs(value) > LARGEST_WRITE_VALUE ? sign * LARGEST_WRITE_VALUE : value;

			file << value << " ";

			value = 0.5 * (v[I][j+1] + v[I][j]);
			sign = value >= 0 ? 1 : -1;
			value = fabs(value) < SMALLEST_WRITE_VALUE ? 0 : value;
			value = fabs(value) > LARGEST_WRITE_VALUE ? sign * LARGEST_WRITE_VALUE : value;

			file << value << " ";

			file << " 0.0\t";
		}
		file << "\n";
	}

	file << "\n";

	file << "SCALARS vorticity float 1" << std::endl;
	file << "LOOKUP_TABLE default" << std::endl;
	for (J = JSTART; J <= JEND; J++)
	{
		j = J;
		for (I = ISTART; I <= IEND; I++)
		{
			i = I;

			value = (1 / dx) * (u[i + 1][J] - u[i][J]) - (1 / dy) * (v[I][j + 1] - v[I][j]);
			sign = value >= 0 ? 1 : -1;
			value = fabs(value) < SMALLEST_WRITE_VALUE ? 0 : value;
			value = fabs(value) > LARGEST_WRITE_VALUE ? sign * LARGEST_WRITE_VALUE : value;

			file << value << "\t";
		}
		file << "\n";
	}

	file << "\n";

	file << "SCALARS streamlines float 1" << std::endl;
	file << "LOOKUP_TABLE default" << std::endl;
	for (J = JSTART; J <= JEND; J++)
	{
		j = J;
		for (I = ISTART; I <= IEND; I++)
		{
			i = I;

			value = 0.5 * dy * (u[i + 1][J] + u[i][J]) - 0.5 * dx * (v[I][j + 1] + v[I][j]);
			sign = value >= 0 ? 1 : -1;
			value = fabs(value) < SMALLEST_WRITE_VALUE ? 0 : value;
			value = fabs(value) > LARGEST_WRITE_VALUE ? sign * LARGEST_WRITE_VALUE : value;

			file << value << "\t";
		}
		file << "\n";
	}

	file << "\n";

	file << "SCALARS alpha float 1" << std::endl;
	file << "LOOKUP_TABLE default" << std::endl;
	for (J = JSTART; J <= JEND; J++)
	{
		j = J;
		for (I = ISTART; I <= IEND; I++)
		{
			i = I;

			value = alpha[I][J];
			sign = value >= 0 ? 1 : -1;
			value = fabs(value) < SMALLEST_WRITE_VALUE ? 0 : value;
			value = fabs(value) > LARGEST_WRITE_VALUE ? sign * LARGEST_WRITE_VALUE : value;

			file << value << "\t";
		}
		file << "\n";
	}

	file << "\n";

	file << "SCALARS mu float 1" << std::endl;
	file << "LOOKUP_TABLE default" << std::endl;
	for (J = JSTART; J <= JEND; J++)
	{
		j = J;
		for (I = ISTART; I <= IEND; I++)
		{
			i = I;

			value = mu[I][J];
			sign = value >= 0 ? 1 : -1;
			value = fabs(value) < SMALLEST_WRITE_VALUE ? 0 : value;
			value = fabs(value) > LARGEST_WRITE_VALUE ? sign * LARGEST_WRITE_VALUE : value;

			file << value << "\t";
		}
		file << "\n";
	}

	file << "\n";

	file << "SCALARS rho float 1" << std::endl;
	file << "LOOKUP_TABLE default" << std::endl;
	for (J = JSTART; J <= JEND; J++)
	{
		j = J;
		for (I = ISTART; I <= IEND; I++)
		{
			i = I;

			value = rho[I][J];
			sign = value >= 0 ? 1 : -1;
			value = fabs(value) < SMALLEST_WRITE_VALUE ? 0 : value;
			value = fabs(value) > LARGEST_WRITE_VALUE ? sign * LARGEST_WRITE_VALUE : value;

			file << value << "\t";
		}
		file << "\n";
	}

	file.close();
}