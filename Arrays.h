#pragma once

#include <iostream>

typedef double* Matrix1D;
typedef double** Matrix2D;

Matrix1D double_1D(int nrows)
{
	double* v = new double[nrows];
	
	for (int i = 0; i < nrows; i++)
	{
		v[i] = 0.0;
	}

	return v;
}

Matrix2D double_2D(int nrows, int ncols)
{
	int nelm = ncols * nrows;

	double** v = new double* [nrows];
	v[0] = new double[nelm];

	for (int i = 1; i < nrows; i++)
	{
		v[i] = v[i - 1] + ncols;
	}

	for (int i = 0; i < nrows; i++)
	{
		for (int j = 0; j < ncols; j++)
		{
			v[i][j] = 0.0;
		}
	}

	return v;
}

void delete_Matrix1D(Matrix1D &m)
{
	delete[](m);
	m = nullptr;
}

void delete_Matrix2D(Matrix2D& m)
{
	delete[](m[0]);
	delete[](m);
	m = nullptr;
}