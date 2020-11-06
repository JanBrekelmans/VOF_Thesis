#pragma once

struct Constants
{
	int NPI = 50;
	int NPJ = 50;
	int EX = 2;

	double XMAX = 1.0;
	double YMAX = 1.0;

	double UNORTH = 0.0;
	double USOUTH = 0.0;
	double VEAST = 0.0;
	double VWEST = 0.0;

	double LARGE = 1e30;
	double SMALL = 1e-30;

	double DT = 5e-4;
	double TOTAL_TIME = 1;

	double OMEGA = 1.5;
	int MAX_IT = 1000;
	double TOLERANCE = 1e-6;

	double GX = 0.0;
	double GY = -100.0;
	
	double MU1 = 0.05;
	double MU2 = 0.01;

	double RHO1 = 2.0;
	double RHO2 = 1.0;

	double SIGMA = 0.0001;

	double COMPRESSION = 1;

	double SMALLEST_WRITE_VALUE = 1e-12;
	double LARGEST_WRITE_VALUE = 1e12;
};