//-------------------------------------------------------------------------------------------------
//					Function in class Output
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
//					Initialize
//-------------------------------------------------------------------------------------------------
Output::Output(double Time, Initialization &Ini)
{
	//get domain dimension from Initialization class
	l = Ini.l; d = Ini.d; width = Ini.width; dx = Ini.dx; dy = Ini.dy; dz = Ini.dz;

}
//-------------------------------------------------------------------------------------------------
//					Output fluid states on grid points
//-------------------------------------------------------------------------------------------------
void Output::OutStates(BaseSolver *CurrentSolver, double Time, Initialization &Ini)
{
	char file_name[50], file_list[50];
	
	//produce output file name
	double Itime = Time*1.0e6;
	strcpy(file_name,"./outdata/time");
	sprintf(file_list, "%d", (int)Itime);
	strcat(file_name, file_list);
	strcat(file_name, ".dat");

	ofstream out(file_name);
	//defining header for tecplot(plot software)
	out<<"title='View'"<<"\n";
	out<<"variables=x, y, z, p, rho, u, v, w, Hybrid_factor"<<"\n";
	out<<"zone t='filed'  i= "<<l+1<<"  j= "<<d+1<<"  k= "<<width+1<<"\n";

	double ttlu = 0.0;
	double rho, p, u, v, w, Hybrid_factor;
	for(int k=4; k<=width+4; k++) 
		for(int j=4; j<=d+4; j++) 
			for(int i=4; i<=l+4; i++) {
				CurrentSolver->OutputStates(i, j, k, rho, p, u, v, w, Hybrid_factor);
				ttlu += u;
				out<<dx*(i-4)+0.5*dx<<" ";
				out<<dy*(j-4)+0.5*dy<<" ";
				out<<dz*(k-4)+0.5*dz<<" ";
				out<<p<<" ";
				out<<rho<<" ";
				out<<u<<" ";
				out<<v<<" ";
				out<<w<<" ";
				out<<Hybrid_factor<<" ";
				out<<"\n";
			}
		out.close();

		cout<<"symmetric error = '"<< ttlu <<"\n";

}
//output error for case 21 problem 8
void  Output::OutError218(BaseSolver *CurrentSolver, double Time, Initialization &Ini)
{
	//produce output file name
	char file_name[50];
	strcpy(file_name,"./outdata/error218.dat");

	double L_1 = 0.0, L_inf = 0.0, ttl_rho = 0.0;

	double beta = 10.0828;
	double beta2 = beta*beta;
	double rho, p, u, v, w, Hybrid_factor;
	for(int k=4; k<=width+4; k++) 
		for(int j=4; j<=d+4; j++) 
			for(int i=4; i<=l+4; i++) {
				double x = (i-4)*dx + 0.5*dx; 
				double y = (j-4)*dy + 0.5*dy; 
				double _x = (x - 5.01);
				double _y = (y - 5.01);
				double r2 = (_x*_x + _y*_y);
				double inter1 = beta2*exp(1.0 - r2)/28.0/pi/pi;
				double rho_theory = pow(1.0 - inter1, 2.5); 
				CurrentSolver->OutputStates(i, j, k, rho, p, u, v, w, Hybrid_factor);
				double error = fabs(rho_theory - rho);
				L_1 += error;
				L_inf = AMAX1(L_inf, error); 
				ttl_rho += rho_theory;
			}

	L_1 = L_1/ttl_rho;

	ofstream out(file_name);
	out<<"L_1 error = "<<L_1<<"  L_inf = "<<L_inf<<"\n";

	out.close();
}
