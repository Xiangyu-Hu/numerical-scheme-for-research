//-------------------------------------------------------------------------------------------------
//							Function in class Initialization
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
Initialization::Initialization()
{
	char Key_word[100];

	ifstream fin("input.dat"); //read input file

	while(!fin.eof()) {
		
		fin>>Key_word;

		//choose the case number
		if(!strcmp(Key_word, "*case_numbers")){
			fin>>case_number>>problem_number;
		}
	
		//resolution in x direction
		if(!strcmp(Key_word, "*resolution_in_x_direction")){
			fin>>l;
		}
		
		//restart file input parameters: 1 read restart file, 0 not
		if(!strcmp(Key_word, "*restart_file")){
	   	 	fin>>restart_file;
		}
	}
	fin.close();

	//set defaults
	length = 1.0; 
	lc = 1.0; //characteristic length
	start_time = 0.0;

	//cases
	if(case_number == 1) { 
		length = 1.0;
		if(problem_number == 3) {
			length = 10.0;
			lc = 10.0;
		}
		if(problem_number == 7) length = 9.0;
		if(problem_number == 9) length = 4.0 + 0.005;
		ld = 0; lw = 0;
		dydx = 1.0; dzdx = 1.0;

		cout<<"\n 1D case with inviscid single fluid of ideal gas. \n";
	}

	//cases
	if(case_number == 21) { 
		length = 1.0;
		if(problem_number == 1 || problem_number == 8) length = 10.0;
		if(problem_number == 5) length = 30.0;
		if(problem_number == 6) length = 1.1;
		ld = 1; lw = 0;
		dydx = 1.0; dzdx = 1.0;

		cout<<"\nTwo-diemnsional cases 1 \n";
	}

	//cases
	if(case_number == 22) { 
		length = 1.0;
		if(problem_number == 1) length = 4.0;
		ld = 4; lw = 0;
		dydx = 1.0; dzdx = 1.0;

		cout<<"\nTwo-diemnsional cases 2 \n";
	}

	//cases
	if(case_number == 23) { 
		length = 1.0;
		if(problem_number == 1) length = 4.0;
		ld = 2; lw = 0;
		dydx = 1.0; dzdx = 1.0;

		cout<<"\nTwo-diemnsional cases 3 \n";
	}

	//cases for testing 2d code
	if(case_number == 24) { 
		length = 1.0;
		if(problem_number == 1) length = 10.0;
		ld = 0; lw = 1;
		dydx = 1.0; dzdx = 1.0;

		cout<<"\nTwo-diemnsional cases 1 \n";
	}

	//cases
	if(case_number == 30) { 
		ld = 1; lw = 1;
		length = 2.0*pi;
		dydx = 1.0; dzdx = 1.0;

		cout<<"\nThree-diemnsional cube with domain size 2*Pi 1 \n";
	}
	//redefine the starting time if restarting form a restart file
	if(restart_file == 1) {
		ifstream fin("singlefluid.rst");
		while(!fin.eof()) {
			fin>>Key_word;
			if(!strcmp(Key_word, "*re_init_time")){
				fin>>start_time;
				break;
			}
		}
		fin.close();
	}

	//domain size
	d = ld == 0 ? 1 : l/ld;
	width = lw == 0 ? 1 : l/lw; 
	cout<<"Domain size in x-direction, length= "<<length<<" \n";
	cout<<"Number spatial steps in x, y and z directions, l= "<<l<<", d = "<<d<<", width = "<<width<<" \n";

	//step-size
	dx = length/double(l);
	dy = dx*dydx;	
	dz = dx*dzdx;	
	
	//adjust l,d for matrix index
	l=l-1; d=d-1, width = width-1;
}
void Initialization::InitialCondition(BaseFluid *ActiveFluid, double rho[Xmax][Ymax][Zmax], double p[Xmax][Ymax][Zmax],
	double u[Xmax][Ymax][Zmax], double v[Xmax][Ymax][Zmax], double w[Xmax][Ymax][Zmax])
{

	//cases
	if(case_number == 1) {

		//problem 0
		//shock tube problem
		if(problem_number == 0) { 
			ActiveFluid->xBl = 2; ActiveFluid->xBr = 2; 
			ActiveFluid->yBu = 0; ActiveFluid->yBd = 0;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						rho[i][j][k] = 1.0; p[i][j][k] = 1.0;
						u[i][j][k] = 1.0; v[i][j][k] = 0.0; w[i][j][k] = 0.0;
						if(i >= int((l+1)*0.4)+4 && i <= int((l+1)*0.6)+4) {
							rho[i][j][k] = 10.0; p[i][j][k] = 1.0;
						}
			}
			end_time = 0.2;
			output_period = 0.2;
		}

		//problem 1
		//shock tube problem
		if(problem_number == 1) { 
			ActiveFluid->xBl = 2; ActiveFluid->xBr = 2; 
			ActiveFluid->yBu = 0; ActiveFluid->yBd = 0;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						rho[i][j][k] = 0.125; p[i][j][k] = 0.1;
						u[i][j][k] = 0.0; v[i][j][k] = 0.0; w[i][j][k] = 0.0;
						if(i <= int((l+1)*0.5)+4) {
							rho[i][j][k] = 1.0; p[i][j][k] = 1.0;
						}
			}
			ActiveFluid->cfl = 0.5;
			end_time = 0.2;
			output_period = 0.2;
		}
		
		//problem 2
		//Lax problem
		if(problem_number == 2) { 
			ActiveFluid->xBl = 2; ActiveFluid->xBr = 2; 
			ActiveFluid->yBu = 0; ActiveFluid->yBd = 0;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			double x_0 = 0.5*length;
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						double x = (i-4)*dx + 0.5*dx; 
						rho[i][j][k] = 0.5; p[i][j][k] = 0.571;
						u[i][j][k] = 0.0; v[i][j][k] = 0.0; w[i][j][k] = 0.0; 
						if(x < x_0) {
							rho[i][j][k] = 0.445; p[i][j][k] = 3.528;
							u[i][j][k] = 0.698; 
						}
			}
			end_time = 0.14;
			output_period = 0.14;
		}

		//problem 3
		//shu-osher problem
		if(problem_number == 3) { 
			ActiveFluid->xBl = 2; ActiveFluid->xBr = 2; 
			ActiveFluid->yBu = 0; ActiveFluid->yBd = 0;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			double x_0 = 0.1*length;
			double x_1 = 0.5*length;
			double wave_number = 50.0 / length;
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						double x = (i-4)*dx + 0.5*dx; 
						rho[i][j][k] = 1.0 + 0.2*sin(wave_number*(x - x_1)); p[i][j][k] = 1.0;
						u[i][j][k] = 0.0; v[i][j][k] = 0.0; w[i][j][k] = 0.0; 
						if(x < x_0) {
							rho[i][j][k] = 3.857143; p[i][j][k] = 10.333333;
							u[i][j][k] = 2.629369;
					}
			}
			end_time = 1.8;
			output_period = 1.8;
		}

		//problem 4
		//two blast wave problem
		if(problem_number == 4) { 
			ActiveFluid->xBl = 0; ActiveFluid->xBr = 0; 
			ActiveFluid->yBu = 0; ActiveFluid->yBd = 0;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			double x_0 = 0.1*length;
			double x_1 = 0.9*length;
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						double x = (i-4)*dx + 0.5*dx; 
						rho[i][j][k] = 1.0; p[i][j][k] = 0.01;
						u[i][j][k] = 0.0; v[i][j][k] = 0.0; w[i][j][k] = 0.0; 
						if(x < x_0) p[i][j][k] = 1000.0;
						if(x > x_1) p[i][j][k] = 100.0;
			}
//			ActiveFluid->CurrentPositivityPreserving = &BaseFluid::PositivityPreserving_1;
			end_time = 0.038;
			output_period = 0.038;
		}

		//problem 5
		//1-2-3 problem
		if(problem_number == 5) { 
			ActiveFluid->xBl = 2; ActiveFluid->xBr = 2; 
			ActiveFluid->yBu = 0; ActiveFluid->yBd = 0;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			double x_0 = 0.5*length;
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						double x = (i-4)*dx + 0.5*dx; 
						rho[i][j][k] = 1.0; p[i][j][k] = 0.4;
						u[i][j][k] = 2.0; v[i][j][k] = 0.0; w[i][j][k] = 0.0; 
						if(x < x_0) u[i][j][k] = -2.0;
			}
//			ActiveFluid->CurrentPositivityPreserving = &BaseFluid::PositivityPreserving_1;
			ActiveFluid->cfl = 0.45;
			end_time = 0.1;
			output_period = 0.1;
		}
		
		//problem 6
		//shock wave reflection problem
		if(problem_number == 6) { 
			ActiveFluid->xBl = 2; ActiveFluid->xBr = 0; 
			ActiveFluid->yBu = 0; ActiveFluid->yBd = 0;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			double x_0 = 0.1*length;
			double x_1 = 0.5*length;
			double wave_number = 50.0 / length;
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						double x = (i-4)*dx + 0.5*dx; 
						rho[i][j][k] = 1.0; p[i][j][k] = 1.0;
						u[i][j][k] = 0.0; v[i][j][k] = 0.0; w[i][j][k] = 0.0; 
						if(x < x_1) {
							rho[i][j][k] = 3.857143; p[i][j][k] = 10.333333;
							u[i][j][k] = 2.629369;
					}
			}
			end_time = 1.0;
			output_period = 0.2;
		}	

		//problem 7
		//LeBlanc shock tube problem
		if(problem_number == 7) { 
			//different specific heat ratio 
			ActiveFluid->Gamma = 5.0/3.0; 
			ActiveFluid->xBl = 0; ActiveFluid->xBr = 0; 
			ActiveFluid->yBu = 0; ActiveFluid->yBd = 0;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			double x_0 = length/3.0;
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						double x = (i-4)*dx + 0.5*dx; 
						rho[i][j][k] = 1.0; p[i][j][k] = 1.0e-1*2.0/3.0;
						u[i][j][k] = 0.0; v[i][j][k] = 0.0; w[i][j][k] = 0.0; 
						if(x > x_0) {
							rho[i][j][k] = 1.0e-3; p[i][j][k] = 1.0e-10*2.0/3.0;
					}
			}
			ActiveFluid->CurrentPositivityPreserving = &BaseFluid::PositivityPreserving_1;
			end_time = 6.0;
			output_period = 6.0;
		}	

		//problem 8
		//1-2-3 problem with vacuum
		if(problem_number == 8) { 
			ActiveFluid->xBl = 2; ActiveFluid->xBr = 2; 
			ActiveFluid->yBu = 0; ActiveFluid->yBd = 0;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			double x_0 = 0.5*length;
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						double x = (i-4)*dx + 0.5*dx; 
						rho[i][j][k] = 1.0; p[i][j][k] = 0.1;
						u[i][j][k] = 2.0; v[i][j][k] = 0.0; w[i][j][k] = 0.0; 
						if(x < x_0) u[i][j][k] = -2.0;
			}

			ActiveFluid->CurrentPositivityPreserving = &BaseFluid::PositivityPreserving_1;
			ActiveFluid->cfl = 0.5;
			end_time = 0.1;
			output_period = 0.1;
		}

		//problem 9
		//Sedov problem with 801 grid points
		if(problem_number == 9) { 
			ActiveFluid->xBl = 2; ActiveFluid->xBr = 2; 
			ActiveFluid->yBu = 0; ActiveFluid->yBd = 0;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			double x_0 = 0.5*length;
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						double x = (i-4)*dx + 0.5*dx; 
						rho[i][j][k] = 1.0; p[i][j][k] = 4.0e-13;
						u[i][j][k] = 0.0; v[i][j][k] = 0.0; w[i][j][k] = 0.0; 
						if(x > x_0 - 0.6*dx && x < x_0 + 0.6*dx) p[i][j][k] = 2.56e8;
			}
			ActiveFluid->CurrentPositivityPreserving = &BaseFluid::PositivityPreserving_1;
			ActiveFluid->cfl = 0.5;
			end_time = 1.0e-3;
			output_period = 1.0e-4;
		}

		//problem 10
		if(problem_number == 10) { 
			ActiveFluid->xBl = 2; ActiveFluid->xBr = 2; 
			ActiveFluid->yBu = 0; ActiveFluid->yBd = 0;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			double x_0 = 0.5*length;
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						double x = (i-4)*dx + 0.5*dx; 
						rho[i][j][k] = 1.0; p[i][j][k] = 10.0;
						u[i][j][k] = 100.0; v[i][j][k] = 0.0; w[i][j][k] = 0.0; 
						if(x > x_0) {
							rho[i][j][k] = 0.1; p[i][j][k] = 1.0;
							u[i][j][k] = -100.0;
						}
			}
			ActiveFluid->CurrentPositivityPreserving = &BaseFluid::PositivityPreserving_1;
			ActiveFluid->cfl = 0.5;
			end_time = 0.004;
			output_period = 0.004;
		}

		//problem 12:broadband sound waves on 128 grid ponits
		if(problem_number == 12) { 
			ActiveFluid->xBl = 1; ActiveFluid->xBr = 1; 
			ActiveFluid->yBu = 0; ActiveFluid->yBd = 0;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			//spectrum
			int hlf_l = 64;
			double *phase; //random phase shift
			phase = new double[hlf_l];
			//creat random noise 
			ntab = 32;
			iv = new long int[ntab];
			Ranils();
			for(int m=0; m<hlf_l; m++)  phase[m] = Ranuls();
			double p0 = 1.0, rho0 = 1.4, c0 = 1.0, u0 = 0.0, epsilon = 1.0e-3, k0 = 12.0;
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						double x = (i-4)*dx + 0.5*dx; 
						double preturb = 0.0; 
						for(int m=0; m<hlf_l; m++) {
							double mp1 = double(m+1);
							preturb += sqrt(exp(-2.0*mp1*mp1/k0/k0)*mp1*mp1*mp1*mp1/k0/k0/k0/k0)*sin(2.0*pi*mp1*(x + phase[m]));
						}
						p[i][j][k] = p0*(1.0 + epsilon*preturb);
						rho[i][j][k] = rho0*pow(p[i][j][k]/p0, 1.0/1.4); 
						u[i][j][k] = u0 + 5.0*sqrt(1.4*p[i][j][k]/rho[i][j][k])/c0; 
						v[i][j][k] = 0.0; w[i][j][k] = 0.0; 
			}
			ActiveFluid->cfl = 0.2;
			end_time = 1.0/6.0;
			output_period = 1.0/6.0;
		}	

		//problem 13:shock-interface interaction
		if(problem_number == 13) { 
			ActiveFluid->xBl = 2; ActiveFluid->xBr = 2; 
			ActiveFluid->yBu = 0; ActiveFluid->yBd = 0;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			double x_1 = 0.15;
			double x_0 = 0.25;
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						double x = (i-4)*dx + 0.5*dx; 
						double y = (j-4)*dy + 0.5*dy; 
						double phi = x-x_0;
						double vlm = Heaviside(phi, 0.025);
						rho[i][j][k] = 1.0*vlm + 1.0; 
						p[i][j][k] = 1.0;
						u[i][j][k] = 0.0;
						v[i][j][k] = 0.0; 
						w[i][j][k] = 0.0; 
						//shock
						if(x < x_1) {
							rho[i][j][k] = 5.268; 
							p[i][j][k] = 41.83;
							u[i][j][k] = 5.752;
						}
			}
			end_time = 0.15;
			output_period = 0.03;
		}	
}

	//cases
	if(case_number == 21) {
		//problem 1
		//vortex translate problem
		if(problem_number == 1) { 
			ActiveFluid->xBl = 1; ActiveFluid->xBr = 1; 
			ActiveFluid->yBu = 1; ActiveFluid->yBd = 1;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			double beta = 5.0;
			double beta2 = beta*beta;
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						double x = (i-4)*dx + 0.5*dx; 
						double y = (j-4)*dy + 0.5*dy; 
						double _x = (x - 5.0);
						double _y = (y - 5.0);
						double r2 = (_x*_x + _y*_y);
						double inter1 = beta2*exp(1.0 - r2)/28.0/pi/pi;
						double inter2 = 0.5*beta*exp(0.5*(1.0 - r2))/pi;
						rho[i][j][k] = pow(1.0 - inter1, 2.5); 
						p[i][j][k] = pow(rho[i][j][k], 1.4);
						u[i][j][k] = 1.0 - inter2*_y;
						v[i][j][k] = 1.0 + inter2*_x; 
						w[i][j][k] = 0.0; 
			}
			end_time = 10.0;
			output_period = 10.0;
		}
	
		//problem 2
		//inviscid double shear layer
		if(problem_number == 2) { 
			ActiveFluid->xBl = 1; ActiveFluid->xBr = 1; 
			ActiveFluid->yBu = 1; ActiveFluid->yBd = 1;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						double x = (i-4)*dx + 0.5*dx; 
						double y = (j-4)*dy + 0.5*dy; 
						rho[i][j][k] = 0.033333333; 
						p[i][j][k] = 1.0e2*rho[i][j][k]/1.4;
						u[i][j][k] = tanh((y - 0.25)/rho[i][j][k]);
						if(y > 0.5) u[i][j][k] = tanh((0.75 - y)/rho[i][j][k]);
						v[i][j][k] = 0.05*sin(2.0*pi*x); 
						w[i][j][k] = 0.0; 
					}
				end_time = 1.8;
				output_period = 0.4;
		}

		//problem 3
		//inviscid Kelvin-Helmholtz instability
		if(problem_number == 3) { 
			ActiveFluid->xBl = 1; ActiveFluid->xBr = 1; 
			ActiveFluid->yBu = 0; ActiveFluid->yBd = 0;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						double num_k = 2.0*pi;
						double x = (i-4)*dx + 0.5*dx; 
						double y = (j-4)*dy + 0.5*dy; 
						rho[i][j][k] = 1.0; 
						p[i][j][k] = 1.0e2*rho[i][j][k]/1.4;
						u[i][j][k] = 0.5 + 0.1*cos(num_k*(x - 0.5))*exp(-num_k*fabs(y - 0.5));
						v[i][j][k] = - 0.1*sin(num_k*(x - 0.5))*exp(-num_k*fabs(y - 0.5));
						if(y < 0.5) u[i][j][k] = - 0.5 - 0.1*cos(num_k*(x - 0.5))*exp(-num_k*fabs(y - 0.5));
						w[i][j][k] = 0.0;
					}
			end_time = 5.0;
			output_period = 1.0;
		}
		//problem 4
		//2D Taylor-Green flow
		if(problem_number == 4) { 
			ActiveFluid->xBl = 1; ActiveFluid->xBr = 1; 
			ActiveFluid->yBu = 1; ActiveFluid->yBd = 1;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						double x = (i-4)*dx + 0.5*dx; 
						double y = (j-4)*dy + 0.5*dy; 
						rho[i][j][k] = 1.0; 
						u[i][j][k] = - cos(2.0*pi*x)*sin(2.0*pi*y);
						v[i][j][k] = sin(2.0*pi*x)*cos(2.0*pi*y); 
						w[i][j][k] = 0.0; 
						p[i][j][k] = 1.0e2/1.4;
					}
			ActiveFluid->mu1 = 0.01; //Reynolds number 100
			ActiveFluid->CurrentViscousHeat_dt = &SingleFluid::ViscousHeat_dt_1;
			ActiveFluid->CurrentPhysicalDissipation = &SingleFluid::PhysicalDissipation_1;

			end_time = 10.0;
			output_period = 5.0;
		}
		//problem 5
		//shock-vortex interaction problem
		if(problem_number == 5) { 
			ActiveFluid->xBl = 2; ActiveFluid->xBr = 3; 
			ActiveFluid->yBu = 3; ActiveFluid->yBd = 3;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			double beta = 5.0;
			double beta2 = beta*beta;
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
					double x = (i-4)*dx + 0.5*dx; 
					double y = (j-4)*dy + 0.5*dy; 
					double _x = (x - 15.0);
					double _y = (y - 15.0);
					double r2 = (_x*_x + _y*_y);
					double Mv = 0.39;
					double rv = 1.0;
					double inter1 = Mv*exp(0.5*(1.0 - r2/rv/rv))/rv;
					double inter2 = 1.0 - 0.2*Mv*Mv*exp(1.0 - r2/rv/rv);
					rho[i][j][k] = pow(inter2, 2.5); 
					u[i][j][k] = -inter1*_y;
					v[i][j][k] = inter1*_x; 
					p[i][j][k] = pow(rho[i][j][k], 1.4)/1.4;
					w[i][j][k] = 0.0; 
					double x_s = 5.0;
					if(x <= x_s) {
						rho[i][j][k] = 1.49826683; 
						u[i][j][k] = 0.429;
						v[i][j][k] = 0.0; 
						p[i][j][k] = 1.77478333/1.4;
					}
				}
			end_time = 16.2;
			output_period = 4.0/1.29;
		}

		//problem 6
		//2d Sedov problem
		if(problem_number == 6) { 
			ActiveFluid->xBl = 0; ActiveFluid->xBr = 3; 
			ActiveFluid->yBu = 3; ActiveFluid->yBd = 0;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						double x = (i-4)*dx + 0.5*dx; 
						double y = (j-4)*dy + 0.5*dy; 
						rho[i][j][k] = 1.0; 
						p[i][j][k] = 4.0e-13;
						u[i][j][k] = 0.0;
						v[i][j][k] = 0.0;
						w[i][j][k] = 0.0;
						if(i==4 && j==4) p[i][j][k] = 0.4*0.244816e6/dx/dy;
					}
			ActiveFluid->cfl = 0.5;
			end_time = 1.0e-3;
			output_period = 1.0e-4;
		}

		//problem 7
		//shock diffraction problem
		if(problem_number == 7) { 
			ActiveFluid->xBl = 217; ActiveFluid->xBr = 3; 
			ActiveFluid->yBu = 3; ActiveFluid->yBd = 3;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			double y_0 = 0.5;
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						double x = (i-4)*dx + 0.5*dx; 
						double y = (j-4)*dy + 0.5*dy; 
						rho[i][j][k] = 1.0; 
						p[i][j][k] = 1.0;
						u[i][j][k] = 0.0;
						v[i][j][k] = 0.0;
						w[i][j][k] = 0.0;
						if(y > y_0 && i <= 4) {
							rho[i][j][k] = 5.268; 
							p[i][j][k] = 41.83;
							u[i][j][k] = 5.752;
					}
					}
			ActiveFluid->cfl = 0.4;
			end_time = 1.0e-1;
			output_period = 2.0e-2;
		}

		//problem 8
		//strong vortex translate problem
		if(problem_number == 8) { 
			ActiveFluid->xBl = 1; ActiveFluid->xBr = 1; 
			ActiveFluid->yBu = 1; ActiveFluid->yBd = 1;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			double beta = 10.0828;
			double beta2 = beta*beta;
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						double x = (i-4)*dx + 0.5*dx; 
						double y = (j-4)*dy + 0.5*dy; 
						double _x = (x - 5.0);
						double _y = (y - 5.0);
						double r2 = (_x*_x + _y*_y);
						double inter1 = beta2*exp(1.0 - r2)/28.0/pi/pi;
						double inter2 = 0.5*beta*exp(0.5*(1.0 - r2))/pi;
						rho[i][j][k] = pow(1.0 - inter1, 2.5); 
						p[i][j][k] = pow(rho[i][j][k], 1.4);
						u[i][j][k] = 1.0 - inter2*_y;
						v[i][j][k] = 1.0 + inter2*_x; 
						w[i][j][k] = 0.0; 
			}
			ActiveFluid->cfl = pow(dx, 5.0/3.0);
			end_time = 1.0e-2;
			output_period = 1.0e-2;
		}

		//problem 9
		//2d Riemann problem #5
		if(problem_number == 9) { 
			ActiveFluid->xBl = 3; ActiveFluid->xBr = 3; 
			ActiveFluid->yBu = 3; ActiveFluid->yBd = 3;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						double x = (i-4)*dx + 0.5*dx; 
						double y = (j-4)*dy + 0.5*dy; 
						double x_c = 0.5;
						double y_c = 0.5;
						rho[i][j][k] = 2.0; 
						p[i][j][k] = 1.0;
						u[i][j][k] = -0.75;
						v[i][j][k] = 0.5;
						if(x > x_c) {
							rho[i][j][k] = 1.0; 
							v[i][j][k] = -0.5;
						}
						if(y < y_c) {
							rho[i][j][k] = 3.0; 
							u[i][j][k] = 0.75;
							if(x < x_c) {
								rho[i][j][k] = 1.0; 
								u[i][j][k] = 0.75;
							}
						}
						w[i][j][k] = 0.0; 
			}
			end_time = 0.23;
			output_period = 0.046;
		}
	}

	//cases
	if(case_number == 22) {
		//problem 1
		//double Mach reflection
		if(problem_number == 1) { 
			ActiveFluid->xBl = 2; ActiveFluid->xBr = 3; 
			ActiveFluid->yBu = 221; ActiveFluid->yBd = 221;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			double x_1 = 1.0/6.0;
			double phi = pi/3.0;
			double slope = tan(phi);
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						rho[i][j][k] = 1.4; 
						p[i][j][k] = 1.0;
						u[i][j][k] = 0.0;
						v[i][j][k] = 0.0; 
						w[i][j][k] = 0.0; 
						double x = (i-4)*dx + 0.5*dx; 
						double y = (j-4)*dy + 0.5*dy; 
						if(y > slope*(x - x_1)) {
							rho[i][j][k] = 8; 
							p[i][j][k] = 140.2/1.2;
							u[i][j][k] = 8.25*sin(phi);
							v[i][j][k] = - 8.25*cos(phi);
						} 			
			}
			ActiveFluid->CurrentPositivityPreserving = &BaseFluid::PositivityPreserving_1;
			end_time = 0.2;
			output_period = 0.05;
		}

		//problem 2
		//Mach 2000 jet
		if(problem_number == 2) { 
			//different specific heat ratio 
			ActiveFluid->Gamma = 5.0/3.0; 
			ActiveFluid->xBl = 2; ActiveFluid->xBr = 3; 
			ActiveFluid->yBu = 3; ActiveFluid->yBd = 0;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			double y_0 = 0.05*length;
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						rho[i][j][k] = 0.5; 
						p[i][j][k] = 0.4127;
						u[i][j][k] = 0.0;
						v[i][j][k] = 0.0; 
						w[i][j][k] = 0.0; 
						double y = (j-4)*dy + 0.5*dy; 
						//in flow
						if(y <= y_0 && i <= 4) {
							rho[i][j][k] = 5.0; 
							u[i][j][k] = 800.0;
						}
			}
			ActiveFluid->cfl = 0.25;
			end_time = 1.0e-3;
			output_period = 1.0e-4;
		}	

	}
	//cases
	if(case_number == 23) {
		//problem 1
		//2d shu-osher problem
		if(problem_number == 1) { 
			ActiveFluid->xBl = 2; ActiveFluid->xBr = 3; 
			ActiveFluid->yBu = 1; ActiveFluid->yBd = 1;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			double x_1 = 0.5;
			double phi = pi/6.0;
			double wnum = 2.0*pi;
			double cs1 = sqrt(1.4);
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						rho[i][j][k] = 5.565; 
						p[i][j][k] = 74.5;
						u[i][j][k] = 7.765;
						v[i][j][k] = 0.0; 
						w[i][j][k] = 0.0; 
						double x = (i-4)*dx + 0.5*dx; 
						double y = (j-4)*dy + 0.5*dy; 
						if(x >= x_1) {
							rho[i][j][k] = 1.0; 
							p[i][j][k] = 1.0;
							u[i][j][k] = - cs1*sin(phi)*cos((x - 1.5)*wnum*cos(phi) 
										 + (y - 1.0)*wnum*sin(phi));
							v[i][j][k] =   cs1*cos(phi)*cos((x - 1.5)*wnum*cos(phi) 
										 + (y - 1.0)*wnum*sin(phi));
						} 
					
			}
			end_time = 0.2;
			output_period = 0.2;
		}

		//problem 2
		//2d shock instabilty problem
		//grid number in x direction 50
		if(problem_number == 2) { 
			ActiveFluid->xBl = 2; ActiveFluid->xBr = 232; 
			ActiveFluid->yBu = 1; ActiveFluid->yBd = 1;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			double x_1 = 6.0/25.0;
			double x_s = 0.5;
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						rho[i][j][k] = 1.0; 
						p[i][j][k] = 1.0;
						u[i][j][k] = 7.101264;
						v[i][j][k] = 0.0; 
						w[i][j][k] = 0.0; 
						double x = (i-4)*dx + 0.5*dx; 
						double y = (j-4)*dy + 0.5*dy; 
						if(x >= x_1) {
							rho[i][j][k] = 5.268; 
							p[i][j][k] = 41.83;
							u[i][j][k] = 1.348;
							if(x < x_1 + 1.5*dx) {
								rho[i][j][k] = x_s + (1.0 - x_s)*5.628; 
								u[i][j][k] = 7.101264/rho[i][j][k];
								p[i][j][k] = (x_s*(3.5 + 0.5*7.101264*7.101264)
											+ (1.0 - x_s)*(3.5*41.83 + 0.5*5.268*1.348*1.348)
											- 0.5*rho[i][j][k]*u[i][j][k]*u[i][j][k])*0.4;
							}
						} 
			}
			ActiveFluid->cfl = 0.5;
			end_time = 20.0;
			output_period = 0.5;
		}	
		//problem 3
		//2d shock drop interaction
		if(problem_number == 3) { 
			ActiveFluid->xBl = 2; ActiveFluid->xBr = 3; 
			ActiveFluid->yBu = 0; ActiveFluid->yBd = 0;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			double x_1 = 0.05;
			double x_0 = 0.25;
			double r_0 = 0.15;
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						double x = (i-4)*dx + 0.5*dx; 
						double y = (j-4)*dy + 0.5*dy; 
						double phi = sqrt((x-x_0)*(x-x_0) + y*y) - r_0;
						double vlm = Heaviside(phi, 0.025);
						rho[i][j][k] = vlm + (1.0 - vlm)*10.0; 
						p[i][j][k] = 1.0;
						u[i][j][k] = -1.0;
						v[i][j][k] = 0.0; 
						w[i][j][k] = 0.0; 
						//shock
						if(x < x_1) {
							rho[i][j][k] = 5.268; 
							p[i][j][k] = 41.83;
							u[i][j][k] = 4.752;
						}
			}
			end_time = 0.15;
			output_period = 0.03;
		}	
		//problem 4
		//2d shock bubble interaction
		if(problem_number == 4) { 
			ActiveFluid->xBl = 2; ActiveFluid->xBr = 3; 
			ActiveFluid->yBu = 0; ActiveFluid->yBd = 0;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			double x_1 = 0.05;
			double x_0 = 0.25;
			double r_0 = 0.15;
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						rho[i][j][k] = 1.0; 
						p[i][j][k] = 1.0;
						u[i][j][k] = -3.0;
						v[i][j][k] = 0.0; 
						w[i][j][k] = 0.0; 
						double x = (i-4)*dx + 0.5*dx; 
						double y = (j-4)*dy + 0.5*dy; 
						//shock
						if(x < x_1) {
							rho[i][j][k] = 5.268; 
							p[i][j][k] = 41.83;
							u[i][j][k] = 2.752;
						}
						//bubble
						if(sqrt((x-x_0)*(x-x_0) + y*y) < r_0) {
							rho[i][j][k] = 0.138; 
						}
			}
			ActiveFluid->CurrentPositivityPreserving = &BaseFluid::PositivityPreserving_1;
			end_time = 0.15;
			output_period = 0.03;
		}	
		
		//problem 5
		//2d viscous shock tube
		if(problem_number == 5) { 
			ActiveFluid->xBl = 4; ActiveFluid->xBr = 4; 
			ActiveFluid->yBu = 0; ActiveFluid->yBd = 4;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			double x_0 = 0.5;
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						rho[i][j][k] = 1.2; 
						p[i][j][k] = 1.2/ActiveFluid->Gamma;
						u[i][j][k] = 0.0;
						v[i][j][k] = 0.0; 
						w[i][j][k] = 0.0; 
						double x = (i-4)*dx + 0.5*dx; 
						double y = (j-4)*dy + 0.5*dy; 
						//shock
						if(x < x_0) {
							rho[i][j][k] = 120.0; 
							p[i][j][k] = 120.0/ActiveFluid->Gamma;
							u[i][j][k] = 0.0;
						}

			}

			ActiveFluid->mu1 = 0.005; //viscosity
			ActiveFluid->mu2 = -0.005*2.0/3.0; //2D Stokes fluid
			ActiveFluid->Pr = 0.73; //Prandtl number
			ActiveFluid->rho_alpha = ActiveFluid->mu1*ActiveFluid->Gamma/ActiveFluid->Pr;
			ActiveFluid->CurrentViscousHeat_dt = &SingleFluid::ViscousHeat_dt_1;
			ActiveFluid->CurrentPhysicalDissipation = &SingleFluid::PhysicalDissipation_1;
//			ActiveFluid->CurrentPositivityPreserving = &BaseFluid::PositivityPreserving_1;

			end_time = 1.0;
			output_period = 0.04;
		}		
}
	
	//cases
	if(case_number == 24) {
		//problem 1
		//vortex translate problem
		if(problem_number == 1) { 
			ActiveFluid->xBl = 1; ActiveFluid->xBr = 1; 
			ActiveFluid->yBu = 0; ActiveFluid->yBd = 0;
			ActiveFluid->zBn = 1; ActiveFluid->zBf = 1;

			//give intial value for the states
			double beta = 5.0;
			double beta2 = beta*beta;
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						double x = (i-4)*dx + 0.5*dx; 
						double z = (k-4)*dz + 0.5*dz; 
						double _x = (x - 5.0);
						double _z = (z - 5.0);
						double r2 = (_x*_x + _z*_z);
						double inter1 = beta2*exp(1.0 - r2)/28.0/pi/pi;
						double inter2 = 0.5*beta*exp(0.5*(1.0 - r2))/pi;
						rho[i][j][k] = pow(1.0 - inter1, 2.5); 
						p[i][j][k] = pow(rho[i][j][k], 1.4);
						u[i][j][k] = 1.0 - inter2*_z;
						v[i][j][k] = 0.0;
						w[i][j][k] = 1.0 + inter2*_x; 
						 
			}
			end_time = 10.0;
			output_period = 10.0;
		}
	
		//problem 2
		//inviscid double shear layer
		if(problem_number == 2) { 
			ActiveFluid->xBl = 1; ActiveFluid->xBr = 1; 
			ActiveFluid->yBu = 1; ActiveFluid->yBd = 1;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						double x = (i-4)*dx + 0.5*dx; 
						double y = (j-4)*dy + 0.5*dy; 
						rho[i][j][k] = 0.033333333; 
						p[i][j][k] = 1.0e2*rho[i][j][k]/1.4;
						u[i][j][k] = tanh((y - 0.25)/rho[i][j][k]);
						if(y > 0.5) u[i][j][k] = tanh((0.75 - y)/rho[i][j][k]);
						v[i][j][k] = 0.05*sin(2.0*pi*x); 
						w[i][j][k] = 0.0; 
					}
				end_time = 1.8;
				output_period = 0.4;
		}

		//problem 3
		//inviscid Kelvin-Helmholtz instability
		if(problem_number == 3) { 
			ActiveFluid->xBl = 1; ActiveFluid->xBr = 1; 
			ActiveFluid->yBu = 0; ActiveFluid->yBd = 0;
			ActiveFluid->zBn = 0; ActiveFluid->zBf = 0;

			//give intial value for the states
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						double num_k = 2.0*pi;
						double x = (i-4)*dx + 0.5*dx; 
						double y = (j-4)*dy + 0.5*dy; 
						rho[i][j][k] = 1.0; 
						p[i][j][k] = 1.0e2*rho[i][j][k]/1.4;
						u[i][j][k] = 0.5 + 0.1*cos(num_k*(x - 0.5))*exp(-num_k*fabs(y - 0.5));
						v[i][j][k] = - 0.1*sin(num_k*(x - 0.5))*exp(-num_k*fabs(y - 0.5));
						if(y < 0.5) u[i][j][k] = - 0.5 - 0.1*cos(num_k*(x - 0.5))*exp(-num_k*fabs(y - 0.5));
						w[i][j][k] = 0.0;
					}
			end_time = 5.0;
			output_period = 1.0;
		}
		//problem 4
		//2D Taylor-Green flow
		if(problem_number == 4) { 
			ActiveFluid->xBl = 1; ActiveFluid->xBr = 1; 
			ActiveFluid->yBu = 0; ActiveFluid->yBd = 0;
			ActiveFluid->zBn = 1; ActiveFluid->zBf = 1;

			//give intial value for the states
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						double x = (i-4)*dx + 0.5*dx; 
						double z = (k-4)*dz + 0.5*dz; 
						rho[i][j][k] = 1.0; 
						u[i][j][k] = - cos(2.0*pi*x)*sin(2.0*pi*z);
						v[i][j][k] = 0.0;
						w[i][j][k] = sin(2.0*pi*x)*cos(2.0*pi*z); 
						p[i][j][k] = 1.0e2/1.4;
					}
			ActiveFluid->mu1 = 0.01; //Reynolds number 100
			ActiveFluid->CurrentViscousHeat_dt = &SingleFluid::ViscousHeat_dt_1;
			ActiveFluid->CurrentPhysicalDissipation = &SingleFluid::PhysicalDissipation_1;

			end_time = 15.0;
			output_period = 1.0;
		}

	}
	
	//cases
	if(case_number == 30) {
		//problem 1
		//inviscid Taylor-Green vortex, weakly compressible flow
		if(problem_number == 1) { 
			ActiveFluid->xBl = 1; ActiveFluid->xBr = 1; 
			ActiveFluid->yBu = 1; ActiveFluid->yBd = 1;
			ActiveFluid->zBn = 1; ActiveFluid->zBf = 1;

			//give intial value for the states
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						double x = (i-4)*dx + 0.5*dx; 
						double y = (j-4)*dy + 0.5*dy; 
						double z = (k-4)*dz + 0.5*dz; 
						
						rho[i][j][k] = 1.0; 
						p[i][j][k] = 1.0e2 + ((cos(2.0*z) + 2)*(cos(2.0*x) + cos(2.0*y)) - 2.0)/16.0;
						u[i][j][k] =   sin(x)*cos(y)*cos(z);
						v[i][j][k] = - cos(x)*sin(y)*cos(z);
						w[i][j][k] = 0.0;
				}
			end_time = 500.0;
			output_period = 1.0;
		}
		//problem 2
		//viscous Taylor-Green vortex, weakly compressible flow
		if(problem_number == 2) { 
			ActiveFluid->xBl = 1; ActiveFluid->xBr = 1; 
			ActiveFluid->yBu = 1; ActiveFluid->yBd = 1;
			ActiveFluid->zBn = 1; ActiveFluid->zBf = 1;

			//give intial value for the states
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						double x = (i-4)*dx + 0.5*dx; 
						double y = (j-4)*dy + 0.5*dy; 
						double z = (k-4)*dz + 0.5*dz; 
						
						rho[i][j][k] = 1.0; 
						p[i][j][k] = 1.0e2 + ((cos(2.0*z) + 2)*(cos(2.0*x) + cos(2.0*y)) - 2.0)/16.0;
						u[i][j][k] =   sin(x)*cos(y)*cos(z);
						v[i][j][k] = - cos(x)*sin(y)*cos(z);
						w[i][j][k] = 0.0;
				}
			//viscous effects
			ActiveFluid->mu1 = 0.01; //Reynolds number 100
			ActiveFluid->CurrentViscousHeat_dt = &SingleFluid::ViscousHeat_dt_1;
			ActiveFluid->CurrentPhysicalDissipation = &SingleFluid::PhysicalDissipation_1;

			end_time = 10.0;
			output_period = 1.0e-4;
			output_period = 5.0;
			}
		//problem 3
		//white noise flow field weakly compressible
		if(problem_number == 3) { 
			ActiveFluid->xBl = 1; ActiveFluid->xBr = 1; 
			ActiveFluid->yBu = 1; ActiveFluid->yBd = 1;
			ActiveFluid->zBn = 1; ActiveFluid->zBf = 1;
			
			//creat random noise 
			ntab = 32;
			iv = new long int[ntab];
			Ranils();

			//give intial value for the states
			for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						double x = (i-4)*dx + 0.5*dx; 
						double y = (j-4)*dy + 0.5*dy; 
						double z = (k-4)*dz + 0.5*dz;

						//three random nubers follows gaussian distribution
						double du, dv, dw;
						Gaussian(du, dv, dw);

						rho[i][j][k] = 1.0; 
						p[i][j][k] = 1.0e2/1.4;
						u[i][j][k] = du;
						v[i][j][k] = dv;
						w[i][j][k] = dw;
				}
			end_time = 100.0;
			output_period = 1.0;
		}
	}
}
//----------------------------------------------------------------------------------------
//			get the random number with uniform distribution in [0, 1]
//----------------------------------------------------------------------------------------
double Initialization::Ranuls()
{
	int j, k;

	//parameters
	long int in1 = 2147483563; 
	long int ik1 = 40014;
	long int iq1 = 53668;
	long int ir1 = 12211;
	long int in2 = 2147483399;
	long int ik2 = 40692;
	long int iq2 = 52774;
	long int ir2 = 3791;
	long int inm1 = in1 - 1;
	long int ndiv = 1 + inm1/ntab;
	double an = 1.0/(double)in1;

//	Linear congruential generactor 1

	k = idum/iq1;
	idum = ik1*(idum - k*iq1) - k*ir1;
	if(idum < 0) idum += in1;

//	Linear congruential generator 2

	k = idum2/iq2;
	idum2 = ik2*(idum2 - k*iq2) - k*ir2;
	if(idum2 < 0) idum2 += in2;

//	Shuffling and subtracting

	j = iy/ndiv;
	iy = iv[j] - idum2;
	iv[j] = idum;

	if(iy < 1) iy += inm1;

	return an*iy;
}
//----------------------------------------------------------------------------------------
//				set the random seed
//----------------------------------------------------------------------------------------
void Initialization::Ranils()
{
	int j, k;

	//for random number generator
	long int in = 2147483563; 
	long int ik = 40014; 
	long int iq = 53668; 
	long int ir = 12211; 

	iseed = rand();

	//	Initial seeds for two random number generators
	idum = iseed + 123456789;
	idum2 = idum;

	//	Load the shuffle table (after 8 warm-ups)
	for( j = ntab+7; j>=0; j--) {
	  k = idum/iq;
	  idum = ik*(idum - k*iq) - k*ir;
	  if(idum < 0) idum = idum + in;
	  if(j < ntab) iv[j] = idum;
	}
	iy = iv[0];
	if(iy < 0) {
		cout<<"iy ="<<iy<<"   Wiener: Random number fails \n"; 
		exit(1);
	}

}
//----------------------------------------------------------------------------------------
//		get three random numbers y1, y2, y3 with guassian distribution 
//		with zero mean and variance one
//		from random numbers uniform distributed in [0, 1]
//----------------------------------------------------------------------------------------
void  Initialization::Gaussian(double &y1, double &y2, double &y3)
{
	double x1, x2, x3, w;
 
	do {
		x1 = 2.0 * Ranuls() - 1.0;
		x2 = 2.0 * Ranuls() - 1.0;
		x3 = 2.0 * Ranuls() - 1.0;
		w = x1 * x1 + x2 * x2 + x3 * x3;
    } while ( w >= 1.0 || w == 0.0);

	w = sqrt( (-2.0 * log( w ) ) / w );
	y1 = x1 * w;
	y2 = x2 * w;
	y3 = x3 * w;
}

