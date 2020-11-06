// VOF.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>

#include "Constants.h"
#include "Grid.h"
#include "State.h"
#include "Solver.h"
#include "Functions.h"
#include "VTK.h"


int main()
{
	Constants constants;
	Grid grid;
	State state;
	Solver_Arrays solver_arrays;

	init_grid(grid, constants);
	init_state(state, grid, constants);
	init_solver_arrays(solver_arrays, constants);

	double time = 0;

	// Set plotting variables and plot initial condition
	int iteration = 0;
	int plot_counter = 1;
	int plot_every_iter = 20;
	write_to_VTK(state, grid, constants, plot_counter++, time);

	while (time < constants.TOTAL_TIME)
	{
		// Update time
		time += constants.DT;
		std::cout << "Time : " << time << "\n";

		// Apply boundary conditions
		bound(state, grid, constants);

		// Calculate normals
		calc_normals(state, grid, constants);

		// Calc ut and vt
		calc_ut(state, grid, constants);
		calc_vt(state, grid, constants);

		calc_alphat(state, grid, constants);
		interface_compression(state, grid, constants);

		// Calculate pressure
		pcoeffs(state, grid, solver_arrays, constants);
		SOR(state.p, solver_arrays, constants);

		// Correct the velocity
		velocity_correct(state, grid, constants);

		// Update variables
		update_states(state, grid, constants);

		// Write the results to VTK
		iteration++;
		if (iteration == plot_every_iter)
		{
			write_to_VTK(state, grid, constants, plot_counter++, time);

			iteration = 0;
		}

		std::cout << "\n";
	}


	// Delete structs
	delete_grid(grid);
	delete_state(state);
	delete_solver_arrays(solver_arrays);
}
