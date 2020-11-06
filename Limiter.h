#pragma once

#include <cmath>

inline double max(double a, double b)
{
	return b > a ? b : a;
}

inline double min(double a, double b)
{
	return b > a ? a : b;
}

inline double van_Leer(double r)
{
	return (r + fabs(r)) / (1 + fabs(r));
}

inline double alpha_limiter(double ap, double an)
{
	//double l = max(pow(1-4*ap*(1-ap)*(1-ap),2),pow(1-4*an*(1-an)*(1-an),2));
	//l = max(1 - l, 0.0);
	//l = min(l, 1.0);
	double L = max(ap * (1 - ap), an * (1 - an));
	L = min(max(L, 0), 1);
	return L;
}