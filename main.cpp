//-------------------------------------------------------------------------------------------------
//		Main program of single compressible flows Xiangyu Hu 2010
//-------------------------------------------------------------------------------------------------
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string.h>
#include <cmath>
#include "glbcls.h"
#include "glbfunc.h"

using namespace std;

int main() {

	double Time;	//computation time
	int N_D_time;	//number of output_period
	BaseSolver *CurrentSolver;


	Initialization Ini; //Geometry and initialization parameters

	//start to compute
	Time = Ini.start_time;
	N_D_time = 0;

	SingleFluidSolver single_fluid_solver(Ini);
	CurrentSolver = &single_fluid_solver;
	single_fluid_solver.InitializeGlobalValues(Time, Ini);

	Output output(Time, Ini);
	output.OutStates(CurrentSolver, Time, Ini);

	while(Time < Ini.end_time) {
		
		// adjust the last output_period
		if(Time + Ini.output_period >=  Ini.end_time) Ini.output_period = Ini.end_time - Time; 
		
		//call the time slover
		CurrentSolver->Time_integeral(Ini, Time, Ini.output_period);
		
		N_D_time = N_D_time + 1;
		
		output.OutStates(CurrentSolver, Time, Ini);
		if(Ini.case_number == 21 && Ini.problem_number == 8) 
			output.OutError218(CurrentSolver, Time, Ini);
	}
	return 0;
}
