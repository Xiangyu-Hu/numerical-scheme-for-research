//-------------------------------------------------------------------------------------------------
//							Function in class BaseSolver
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
//-------------------------------------------------------------------------------------------------
//								Initialize
//-------------------------------------------------------------------------------------------------
BaseSolver::BaseSolver()
{
	Ite = 0;

	char Key_word[100];

	ifstream fin("input.dat"); //read input file

	while(!fin.eof()) {
		
		fin>>Key_word;
		
		//Read domain size and dimenionless parameters
		if(!strcmp(Key_word, "*checking_points")){
			fin>>CheckPoints[0]>>CheckPoints[1];
		
			//output on screen
			cout<<"\nStepwise output parameters\n";
			cout<<"Output conservatives every "<<CheckPoints[0]<<" time steps\n";
			cout<<"Output restart files every "<<CheckPoints[1]<<" time steps\n\n";
		}
	}
	fin.close();
}
//-------------------------------------------------------------------------------------------------
//								Initialize single-fluid solver
//-------------------------------------------------------------------------------------------------
SingleFluidSolver::SingleFluidSolver(Initialization &Ini)
{
	//get domain dimension from Initialization class
	l = Ini.l; d = Ini.d; width = Ini.width; 
	dx = Ini.dx; dy = Ini.dy; dz = Ini.dz;
	CurrentFluid = &fluid;
	fluid.FluidInitializer(Ini);
}
//-------------------------------------------------------------------------------------------------
//				Initialize total entropy before the first step
//-------------------------------------------------------------------------------------------------
void SingleFluidSolver::InitializeGlobalValues(double Time, Initialization &Ini)
{
	double Itime; //current time
	Itime = Time + 1.0e-13; //produce output file name
	
	//total entropy and kinetic energy
	TotalKE0 = 0.0; Max_U0 = 0.0; TotalVolume = 0.0;
	
		//summation
	for(int k=4; k<=width+4; k++) 
		for(int j=4; j<=d+4; j++) 
			for(int i=4; i<=l+4; i++) {
			//Kinetic energy
			double KE = 0.5*(CurrentFluid->U[1][i][j][k]*CurrentFluid->U[1][i][j][k]
			+ CurrentFluid->U[2][i][j][k]*CurrentFluid->U[2][i][j][k]
			+ CurrentFluid->U[3][i][j][k]*CurrentFluid->U[3][i][j][k])/CurrentFluid->U[0][i][j][k];
			TotalKE0 += KE;
			TotalVolume += dx*dy*dz;
			Max_U0 = AMAX1(Max_U0, sqrt(CurrentFluid->u[i][j][k]*CurrentFluid->u[i][j][k] 
									  + CurrentFluid->v[i][j][k]*CurrentFluid->v[i][j][k] 
									  + CurrentFluid->w[i][j][k]*CurrentFluid->w[i][j][k]));
	}
	
	TotalKE0 = TotalKE0*dx*dy*dz/TotalVolume;
	//file for conservation properties
	ofstream out("./outdata/glbvalues.dat"); 
	//header needed for tecplot(plot software)
	out<<"title='View'"<<"\n";
	out<<"variables=time, TotalKE, Max_U"<<"\n";
	out<<Itime<<"  "<<TotalKE0<<"  "<<Max_U0<<"\n";
	out.close();
}
//-------------------------------------------------------------------------------------------------
//				Output total entropy before the first step
//-------------------------------------------------------------------------------------------------
void SingleFluidSolver::Out_GlobalValues(double Time, Initialization &Ini)
{
	double Itime; //current time
	Itime = Time; //produce output file name
	
	//total entropy and kinetic energy
	double TotalKE = 0.0;
	double Max_U = 0.0; 

	//summation
	for(int k=4; k<=width+4; k++) 
		for(int j=4; j<=d+4; j++) 
			for(int i=4; i<=l+4; i++) {
			//Kinetic energy
			double KE = 0.5*(CurrentFluid->U[1][i][j][k]*CurrentFluid->U[1][i][j][k]
			+ CurrentFluid->U[2][i][j][k]*CurrentFluid->U[2][i][j][k]
			+ CurrentFluid->U[3][i][j][k]*CurrentFluid->U[3][i][j][k])/CurrentFluid->U[0][i][j][k];
			TotalKE += KE;
			Max_U = AMAX1(Max_U, sqrt(CurrentFluid->u[i][j][k]*CurrentFluid->u[i][j][k] 
									  + CurrentFluid->v[i][j][k]*CurrentFluid->v[i][j][k] 
									  + CurrentFluid->w[i][j][k]*CurrentFluid->w[i][j][k]));
	}

	//normalized variation
	TotalKE = TotalKE*dx*dy*dz/TotalVolume; 

	//file for conservation properties
	ofstream out("./outdata/glbvalues.dat", ios::out | ios::app); 
	out<<Itime<<"  "<<TotalKE<<"  "<<Max_U<<"\n";
	out.close();

}
//-------------------------------------------------------------------------------------------------
//				Out Conservatives for reinitiation
//-------------------------------------------------------------------------------------------------
void SingleFluidSolver::Out_reinitiation(double Time, Initialization &Ini)
{
	int i,j;
	char file_name[50], file_list[50];
	
	//restart file name
	strcpy(file_name,"singlefluid");
	sprintf(file_list, "%d", Ite/CheckPoints[1]);
	strcat(file_name, file_list);
	strcat(file_name, ".rst");
	ofstream out(file_name);

	//out reinitiation time with dimension
	out<<"*re_init_time"<<"\n";
	out<<Time<<"\n";

	//out fluid1 state and conservatives without dimension
	out<<"*conservatives"<<"\n";
	for(i=0; i<=l+8; i++) 
		for(j=0; j<=d+8; j++) 
			for(int k=0; k<=width+8; k++) 
				for(int n=0; n<5; n++) 
					out<<CurrentFluid->U[n][i][j][k]<<" ";
	out.close();
}
//-------------------------------------------------------------------------------------------------
//								Time integral for single fluid solver
//-------------------------------------------------------------------------------------------------
void SingleFluidSolver::Time_integeral(Initialization &Ini, double &Time, double output_period)
{
	Integeral_time = 0.0;
	while(Integeral_time < output_period) {

		double dt =  CurrentFluid->Get_dt();
		//fitting output period
		if(Integeral_time + dt > output_period) dt = output_period - Integeral_time;
		Integeral_time = Integeral_time + dt;
		Ite = Ite +1;
		Time = Time + dt;

		//screeen output
		if(Ite % CheckPoints[0] == 0) cout<<"N="<<Ite<<" Time: "
							 <<Time<<"	dt: "<<dt<<"\n";

		//solved the fluid with 3rd order Runge-Kutta method
		CurrentFluid->Step_rk3(Time, dt, 1);
		CurrentFluid->Step_rk3(Time, dt, 2);
		CurrentFluid->Step_rk3(Time, dt, 3);

		if(Ite % CheckPoints[0] == 0) Out_GlobalValues(Time, Ini);
		if(Ite % CheckPoints[1] == 0) Out_reinitiation(Time, Ini);
	}
	cout << "Excution time =" << CurrentFluid->exec_time - CurrentFluid->time_start << "\n";

}
//-------------------------------------------------------------------------------------------------
//				Get the state at a grid point
//-------------------------------------------------------------------------------------------------
void SingleFluidSolver::OutputStates(int i, int j, int k, double &rho, double &p, 
									 double &u, double &v, double &w, double &Hybrid_factor)
{
	CurrentFluid->OutputStates(i, j, k, rho, p, u, v, w, Hybrid_factor);
}
