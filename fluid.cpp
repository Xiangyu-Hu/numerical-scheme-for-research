//-----------------------------------------------------------------------------------------------
//								Function in class Fluid
//-------------------------------------------------------------------------------------------------
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string.h>
#include <cmath>
#include <ctime>
#include "glbcls.h"
#include "glbfunc.h"

using namespace std;
//-------------------------------------------------------------------------------------------------
//		BaseFluid constructor
//-------------------------------------------------------------------------------------------------
BaseFluid::BaseFluid()
{
	char Key_word[100];

	ifstream fin("input.dat"); //read input file

	while(!fin.eof()) {
		
		fin>>Key_word;
		//the interpolation method
		if(!strcmp(Key_word, "*interpolation_method")){
	   	 	fin>>interpolation_method;
		}

		//the numerical flux method
		if(!strcmp(Key_word, "*flux_method")){
	   	 	fin>>flux_method;
		}
	}
	fin.close();

	//constant for ideal gas EOS
	Gamma = 1.4; //specific heat ratio 
	_Gamma = 0.4; //Gamma -1
	_Gamma1 = 2.5; //1 / (Gamma -1)
	cfl = 0.5;

	//default values
	g_x = 0.0; g_y = 0.0;
	mu1 = 0.0; mu2 = 0.0; Pr = 1.0; rho_alpha = mu1*Gamma/Pr;
	CurrentViscousHeat_dt = &BaseFluid::ViscousHeat_dt_0;
	CurrentPhysicalDissipation = &BaseFluid::PhysicalDissipation_0;
	CurrentPositivityPreserving = &BaseFluid::PositivityPreserving_0;
}
//-------------------------------------------------------------------------------------------------
//		Boundary condition
//		0 Rigid wall condition
//		1 Peroidic condition
//		2 Inflow condition
//		3 Outflow condition
//-------------------------------------------------------------------------------------------------
void BaseFluid::Boundary(double UI[Emax][Xmax][Ymax][Zmax])
{
	//x direction
	for(int k=4; k<=width+4; k++) 
		for(int j=4; j<=d+4; j++) {
		//for left side
		for(int i=3; i>=0; i--){
			switch(xBl) {
			//the rigid wall conditions 	
			case 0:
				for(int n=0; n<5; n++) 
					UI[n][i][j][k] = UI[n][7-i][j][k];
				UI[1][i][j][k] *= -1.0;
				break;
			//the peroidic conditions	
			case 1:
				for(int n=0; n<5; n++) 
					UI[n][i][j][k] = UI[n][l+i+1][j][k];
				break;
			//the inflow conditions 	
			case 2:
				for(int n=0; n<5; n++) UI[n][i][j][k] = U[n][i][j][k];
				break;
			//the outflow conditions	
			case 3:
				for(int n=0; n<5; n++) 
					UI[n][i][j][k] = UI[n][i+1][j][k];
				break;
			//the no-slip conditions	
			case 4:
				for(int n=0; n<5; n++) 
					UI[n][i][j][k] = UI[n][7-i][j][k];
				UI[1][i][j][k] *= -1.0;
				UI[2][i][j][k] *= -1.0;
				UI[3][i][j][k] *= -1.0;
				break;
			//boundary condition for case 21 problem 7
			case 217:
			double y = (j-4)*dy + 0.5*dy; 
			for(int n=0; n<5; n++) UI[n][i][j][k] = U[n][i][j][k];
			if(y <= 0.5) {
				for(int n=0; n<5; n++) UI[n][i][j][k] = UI[n][7-i][j][k];
				UI[1][i][j][k] = -UI[1][7-i][j][k];
				}
				break;
			}//end-switch
		}
		
		//for right side
		for(int i=l+5; i<=l+8; i++){
			switch(xBr) {
			//the rigid wall conditions 	
			case 0:
				for(int n=0; n<5; n++) 
					UI[n][i][j][k] = UI[n][2*l+9-i][j][k];
				UI[1][i][j][k] *= -1.0;
				break;
			//the peroidic conditions	
			case 1:
				for(int n=0; n<5; n++) 
					UI[n][i][j][k] = UI[n][i-l-1][j][k];
				break;
			//the inflow conditions 	
			case 2:
				for(int n=0; n<5; n++) UI[n][i][j][k] = U[n][i][j][k];
				break;
			//the outflow conditions	
			case 3:
				for(int n=0; n<5; n++) 
					UI[n][i][j][k] = UI[n][i-1][j][k];
				break;
			//the no-slip conditions 	
			case 4:
				for(int n=0; n<5; n++) 
					UI[n][i][j][k] = UI[n][2*l+9-i][j][k];
				UI[1][i][j][k] *= -1.0;
				UI[2][i][j][k] *= -1.0;
				UI[3][i][j][k] *= -1.0;
				break;
			//user defined boundary condition for case 23 problem 2	
			case 232:
				for(int n=0; n<5; n++){ 
					if(n != 1) UI[n][i][j][k] = UI[n][i-1][j][k];
				}
				break;
			}//end-switch
		}
	}

	//y direction
	for(int k=4; k<=width+4; k++)
		for(int i=0; i<=l+8; i++){
		//the bottom side
		for(int j=3; j>=0; j--) {
			switch(yBd) {
			//the rigid wall conditions 	
			case 0:
				for(int n=0; n<5; n++)
					UI[n][i][j][k] = UI[n][i][7-j][k];
				UI[2][i][j][k] *= -1.0;
				break;
			//the peroidic conditions	
			case 1:
				for(int n=0; n<5; n++) 
					UI[n][i][j][k] = UI[n][i][j+d+1][k];
				break;
			//the inflow conditions 	
			case 2:
				for(int n=0; n<5; n++) UI[n][i][j][k] = U[n][i][j][k];
				break;
				//the outflow conditions	
			case 3:
				for(int n=0; n<5; n++) 
					UI[n][i][j][k] = UI[n][i][j+1][k];
				break;
			//the no-slip conditions 	
			case 4:
				for(int n=0; n<5; n++)
					UI[n][i][j][k] = UI[n][i][7-j][k];
				UI[1][i][j][k] *= -1.0;
				UI[2][i][j][k] *= -1.0;
				UI[3][i][j][k] *= -1.0;
				break;
			//user defined boundary condition for case 22 problem 1	
			case 221:
				double phi = pi/3.0;
				double x_1 = 1.0/6.0;
				double x = (i-4)*dx + 0.5*dx; 
				for(int n=0; n<5; n++) 
					UI[n][i][j][k] = UI[n][i][7-j][k];
				UI[2][i][j][k] = -UI[2][i][7-j][k];
				if(x < x_1) {
					rho[i][j][k] = 1.4*120.0/21.0; 
					u[i][j][k] = 99.0/12.0*sin(phi);
					v[i][j][k] = - 99.0/12.0*cos(phi);
					w[i][j][k] = 0.0;
					p[i][j][k] = 140.2/1.2;
					GetU(i, j, k, UI);
				} 
				break;
			}//end-switch
		}
				
			//the top side
			for(int j=d+5; j<=d+8; j++) {
				switch(yBu) {
				//the rigid wall conditions 	
				case 0:
					for(int n=0; n<5; n++) 
						UI[n][i][j][k] = UI[n][i][2*d+9-j][k];
					UI[2][i][j][k] *= -1.0;
					break;
					//the peroidic conditions	
				case 1:
					for(int n=0; n<5; n++) 
						UI[n][i][j][k] = UI[n][i][j-d-1][k];
					break;
				//the inflow conditions 	
				case 2:
					for(int n=0; n<5; n++) UI[n][i][j][k] = U[n][i][j][k];
					break;
				//the outflow conditions	
				case 3:
					for(int n=0; n<5; n++) 
						UI[n][i][j][k] = UI[n][i][j-1][k];
					break;
				//the no-slip conditions 	
				case 4:
					for(int n=0; n<5; n++) 
						UI[n][i][j][k] = UI[n][i][2*d+9-j][k];
					UI[1][i][j][k] *= -1.0;
					UI[2][i][j][k] *= -1.0;
					UI[3][i][j][k] *= -1.0;
					break;
					//user defined boundary condition for case 22 problem 1	
				case 221:
					double phi = pi/3.0;
					double slope = tan(phi);
					double x_1 = 1.0/6.0 + run_time*10.0/sin(phi);
					double x = (i-4)*dx + 0.5*dx; 
					double y = (j-4)*dy + 0.5*dy; 
					rho[i][j][k] = 1.4; 
					p[i][j][k] = 1.0;
					u[i][j][k] = 0.0;
					v[i][j][k] = 0.0; 
					w[i][j][k] = 0.0;
					if(y > slope*(x - x_1)) {
						rho[i][j][k] = 1.4*120.0/21.0; 
						u[i][j][k] = 99.0/12.0*sin(phi);
						v[i][j][k] = - 99.0/12.0*cos(phi);
						p[i][j][k] = 140.2/1.2;
					}
					GetU(i, j, k, UI);
					break;
				}//end-switch
			}
		}
	
	//z direction
	for(int j=0; j<=d+8; j++)
		for(int i=0; i<=l+8; i++){
			//the near side
			for(int k=3; k>=0; k--) {
				switch(zBn) {
				//the rigid wall conditions 	
				case 0:
					for(int n=0; n<5; n++)
						UI[n][i][j][k] = UI[n][i][j][7-k];
					UI[3][i][j][k] *= -1.0;
					break;
				//the peroidic conditions	
				case 1:
					for(int n=0; n<5; n++) 
						UI[n][i][j][k] = UI[n][i][j][k+width+1];
					break;
				//the inflow conditions 	
				case 2:
					for(int n=0; n<5; n++) UI[n][i][j][k] = U[n][i][j][k];
					break;
					//the outflow conditions	
				case 3:
					for(int n=0; n<5; n++) 
						UI[n][i][j][k] = UI[n][i][j][k+1];
					break;
				//the no-slip conditions 	
				case 4:
					for(int n=0; n<5; n++)
						UI[n][i][j][k] = UI[n][i][j][7-k];
					UI[1][i][j][k] *= -1.0;
					UI[2][i][j][k] *= -1.0;
					UI[3][i][j][k] *= -1.0;
					break;
					}//end-switch
			}
				
			//the far side
			for(int k=width+5; k<=width+8; k++) {
				switch(zBf) {
				//the rigid wall conditions 	
				case 0:
					for(int n=0; n<5; n++) 
						UI[n][i][j][k] = UI[n][i][j][2*width+9-k];
					UI[3][i][j][k] *= -1.0;
					break;
					//the peroidic conditions	
				case 1:
					for(int n=0; n<5; n++) 
						UI[n][i][j][k] = UI[n][i][j][k-width-1];
					break;
				//the inflow conditions 	
				case 2:
					for(int n=0; n<5; n++) UI[n][i][j][k] = U[n][i][j][k];
					break;
				//the outflow conditions	
				case 3:
					for(int n=0; n<5; n++) 
						UI[n][i][j][k] = UI[n][i][j][k-1];
					break;
				//the no-slip conditions 	
				case 4:
					for(int n=0; n<5; n++) 
						UI[n][i][j][k] = UI[n][i][j][2*width+9-k];
					UI[1][i][j][k] *= -1.0;
					UI[2][i][j][k] *= -1.0;
					UI[3][i][j][k] *= -1.0;
					break;
				}//end-switch
			}
		}
}
//-------------------------------------------------------------------------------------------------
//		Obtain states from conservative variables
//-------------------------------------------------------------------------------------------------
void BaseFluid::OutputStates(int i, int j, int k, double &denisty, double &pressure, 
							 double &uu, double &vv,  double &ww, double &Hybrid_factor)
{
	denisty	=	rho[i][j][k];
	uu	=	u[i][j][k];
	vv	=	v[i][j][k];
	ww	=	w[i][j][k];
	pressure	=	p[i][j][k];
				
	double EN_x = 0.5*(w[i][j+1][k] - w[i][j-1][k])/dy
				- 0.5*(v[i][j][k+1] - v[i][j][k-1])/dz;
	double EN_y = 0.5*(u[i][j][k+1] - u[i][j][k-1])/dz
				- 0.5*(w[i+1][j][k] - w[i-1][j][k])/dx;
	double EN_z = 0.5*(v[i+1][j][k] - v[i-1][j][k])/dx
				- 0.5*(u[i][j+1][k] - u[i][j-1][k])/dy;
	Hybrid_factor = 1.0 - hybrid_indcator[1][i][j][k]/(hybrid_indcator[0][i][j][k] + 1.0e-15);

}
//-------------------------------------------------------------------------------------------------
//		Obtain flux from cell wall in x direction
//-------------------------------------------------------------------------------------------------
void BaseFluid::CharacteristicF_x(int i, int j, int k, double UI[Emax][Xmax][Ymax][Zmax])
{
	//obtain the eigen vectors
	double eigen_l[Emax][Emax], eigen_r[Emax][Emax], eigen_value[Emax], ave_rho1;
	RoeAverage_x(i, j, k, eigen_l, eigen_r, eigen_value, ave_rho1);
	//Scalar u, f, f_plus, f_minus
	double uf[10], ff[10], pp[10], mm[10], f_flux, _p[Emax][Emax];
	//construct the right value & the left value scalar equations by characteristic reduction			
	// at i+1/2 in x direction
	//prepare for hybrid method
	double du[5], f2[5], du_s[5], f2_s[5];
	//interpolation and computing details
	for(int n=0; n<5; n++){
		for(int m=0; m<stencil_size; m++){
			uf[m] = UI[n][m+i-stencil_P][j][k];
			ff[m] = F_x[n][m+i-stencil_P][j][k];
		}
		du[n] = (*CurrentLinearInterpolation_du)(&uf[stencil_P], dx);
		f2[n] = (*CurrentLinearInterpolation_f2)(&ff[stencil_P], dx);
	}
	//linear decomposition
	for(int n=0; n<5; n++){
		du_s[n] = 0.0; f2_s[n] = 0.0;
		for(int n1=0; n1<5; n1++){
				du_s[n] += du[n1]*eigen_l[n][n1];
				f2_s[n] += f2[n1]*eigen_l[n][n1];
		}
	}
	for(int n=0; n<5; n++){
		//lax-friedrichs
		double eigen_local_max = 0.0; 
		for(int m=0; m<stencil_size; m++) 
			eigen_local_max = AMAX1(eigen_local_max,  fabs(eigen_local[n][m+i-stencil_P][j][k]));
		double ul = eigen_local[n][i][j][k];
		double ur = eigen_local[n][i+1][j][k];
		if (ul < ur && ul*ur < 0.0) eigen_value[n] = 0.5*eigen_value[n] + 0.5*(ur - ul);
		double artificial_viscosity = Roe_type*eigen_value[n]
									+ LLF_type*eigen_local_max + GLF_type*eigen_max[n];
		//hybrid method
		hybrid_indcator[0][i][j][k] += 1.0;
		double sigma = du_s[n]*du_s[n]*ave_rho1*ave_rho1;
		f_flux = f2_s[n] + eigen_value[n] * du_s[n];
		if( sigma > hybrid_epsilon) {
			hybrid_indcator[1][i][j][k] += 1.0;
			double nonlinear_epsilon = n == 3 || n == 4 ? 1.e-6 : 1.0e-3;
			double alpha = Heaviside(sigma - hybrid_epsilon - nonlinear_epsilon, nonlinear_epsilon);
			//			double alpha = 1.0 - Wenland(sigma - hybrid_epsilon, nonlinear_epsilon);

			//characteristic decomposition
			for(int m=0; m<stencil_size; m++){
				uf[m] = 0.0;
				ff[m] = 0.0;
				for(int n1=0; n1<5; n1++){
					uf[m] += UI[n1][m+i-stencil_P][j][k]*eigen_l[n][n1];
					ff[m] += F_x[n1][m+i-stencil_P][j][k]*eigen_l[n][n1];
				}
				pp[m] = 0.5*(ff[m] + artificial_viscosity*uf[m]);  
				mm[m] = 0.5*(ff[m] - artificial_viscosity*uf[m]);  
			}
			// calculate the scalar numerical flux at x direction
			f_flux = (*CurrentInterpolation_P)(&pp[stencil_P], alpha)
				+ (*CurrentInterpolation_M)(&mm[stencil_P], alpha);
		}
		// get Fp
		for(int n1=0; n1<5; n1++)
			_p[n][n1] = f_flux*eigen_r[n1][n];
	}
	// reconstruction the F_x terms
	for(int n=0; n<5; n++){
		F_x_wall[n][i][j][k] = 0.0;
		for(int n1=0; n1<5; n1++) 
			F_x_wall[n][i][j][k] += _p[n1][n];
	}

}
//-------------------------------------------------------------------------------------------------
//		Obtain flux from cell wall in y direction
//-------------------------------------------------------------------------------------------------
void BaseFluid::CharacteristicF_y(int i, int j, int k, double UI[Emax][Xmax][Ymax][Zmax])
{
	//obtain the eigen vectors
	double eigen_l[Emax][Emax], eigen_r[Emax][Emax], eigen_value[Emax], ave_rho1;
	RoeAverage_y(i, j, k, eigen_l, eigen_r, eigen_value, ave_rho1);
	//Scalar u, f, f_plus, f_minus
	double uf[10], ff[10], pp[10], mm[10], f_flux, _p[Emax][Emax];
	//construct the right value & the left value scalar equations by characteristic reduction
	// at j+1/2 in y direction
	//prepare for hybrid method
	double du[5], f2[5], du_s[5], f2_s[5];
	//interpolation and computing details
	for(int n=0; n<5; n++){
		for(int m=0; m<stencil_size; m++){
			uf[m] = UI[n][i][m+j-stencil_P][k];
			ff[m] = F_y[n][i][m+j-stencil_P][k];
		}
		du[n] = (*CurrentLinearInterpolation_du)(&uf[stencil_P], dy);
		f2[n] = (*CurrentLinearInterpolation_f2)(&ff[stencil_P], dy);
	}
	//linear decomposition
	for(int n=0; n<5; n++){
		du_s[n] = 0.0; f2_s[n] = 0.0;
		for(int n1=0; n1<5; n1++){
				du_s[n] += du[n1]*eigen_l[n][n1];
				f2_s[n] += f2[n1]*eigen_l[n][n1];
		}
	}
	for(int n=0; n<5; n++){	
		//lax-friedrichs
		double eigen_local_max = 0.0; 
		for(int m=0; m<stencil_size; m++) 
			eigen_local_max = AMAX1(eigen_local_max,  fabs(eigen_local[n][i][m+j-stencil_P][k]));
		double ul = eigen_local[n][i][j][k];
		double ur = eigen_local[n][i][j+1][k];
		if (ul < ur && ul*ur < 0.0) eigen_value[n] = 0.5*eigen_value[n] + 0.5*fabs(ul - ur);
		double artificial_viscosity = Roe_type*eigen_value[n]
									+ LLF_type*eigen_local_max + GLF_type*eigen_max[n];
		//hybrid method
		hybrid_indcator[0][i][j][k] += 1.0;
		if(du_s[n]*du_s[n]*ave_rho1*ave_rho1 < hybrid_epsilon) {
			f_flux = f2_s[n] + eigen_value[n]*du_s[n]; 
		} else {
			hybrid_indcator[1][i][j][k] += 1.0;
			//characteristic decomposition
			for(int m=0; m<stencil_size; m++){
				uf[m] = 0.0;
				ff[m] = 0.0;
				for(int n1=0; n1<5; n1++){
					uf[m] += UI[n1][i][m+j-stencil_P][k]*eigen_l[n][n1];
					ff[m] += F_y[n1][i][m+j-stencil_P][k]*eigen_l[n][n1];
				}
				pp[m] = 0.5*(ff[m] + artificial_viscosity*uf[m]);  
				mm[m] = 0.5*(ff[m] - artificial_viscosity*uf[m]);  
			}
				
			// calculate the scalar numerical flux at x direction
			f_flux = (*CurrentInterpolation_P)(&pp[stencil_P], ave_rho1)
				+ (*CurrentInterpolation_M)(&mm[stencil_P], ave_rho1);
		}
			// get Gp
		for(int n1=0; n1<5; n1++)
			_p[n][n1] = f_flux*eigen_r[n1][n];
	}
	// reconstruction the F_y flux terms
	for(int n=0; n<5; n++){
		F_y_wall[n][i][j][k] = 0.0;
		for(int n1=0; n1<5; n1++)
			F_y_wall[n][i][j][k] += _p[n1][n];
	}
}
//-------------------------------------------------------------------------------------------------
//		Obtain flux from cell wall in z direction
//-------------------------------------------------------------------------------------------------
void BaseFluid::CharacteristicF_z(int i, int j, int k, double UI[Emax][Xmax][Ymax][Zmax])
{
	//obtain the eigen vectors
	double eigen_l[Emax][Emax], eigen_r[Emax][Emax], eigen_value[Emax], ave_rho1;
	RoeAverage_z(i, j, k, eigen_l, eigen_r, eigen_value, ave_rho1);
	//Scalar u, f, f_plus, f_minus
	double uf[10], ff[10], pp[10], mm[10], f_flux, _p[Emax][Emax];
	//construct the right value & the left value scalar equations by characteristic reduction			
	// at k+1/2 in z direction
	//prepare for hybrid method
	double du[5], f2[5], du_s[5], f2_s[5];
	//interpolation and computing details
	for(int n=0; n<5; n++){
		for(int m=0; m<stencil_size; m++){
			uf[m] = UI[n][i][j][m+k-stencil_P];
			ff[m] = F_z[n][i][j][m+k-stencil_P];
		}
		du[n] = (*CurrentLinearInterpolation_du)(&uf[stencil_P], dz);
		f2[n] = (*CurrentLinearInterpolation_f2)(&ff[stencil_P], dz);
	}
	//linear decomposition
	for(int n=0; n<5; n++){
		du_s[n] = 0.0; f2_s[n] = 0.0;
		for(int n1=0; n1<5; n1++){
				du_s[n] += du[n1]*eigen_l[n][n1];
				f2_s[n] += f2[n1]*eigen_l[n][n1];
		}
	}
	for(int n=0; n<5; n++){
		//lax-friedrichs
		double eigen_local_max = 0.0; 
		for(int m=0; m<stencil_size; m++) 
			eigen_local_max = AMAX1(eigen_local_max,  fabs(eigen_local[n][i][j][m+k-stencil_P]));
		//entropy fix
		double ul = eigen_local[n][i][j][k];
		double ur = eigen_local[n][i][j][k+1];
		if (ul < ur && ul*ur < 0.0) eigen_value[n] = 0.5*eigen_value[n] + 0.5*fabs(ul - ur);
		double artificial_viscosity = Roe_type*eigen_value[n]
									+ LLF_type*eigen_local_max + GLF_type*eigen_max[n];
		//hybrid method
		hybrid_indcator[0][i][j][k] += 1.0;
		if(du_s[n]*du_s[n]*ave_rho1*ave_rho1 < hybrid_epsilon) {
			f_flux = f2_s[n] + eigen_value[n]*du_s[n]; 
		} else {
			hybrid_indcator[1][i][j][k] += 1.0;		
			//characteristic decomposition
			for(int m=0; m<stencil_size; m++){
				uf[m] = 0.0;
				ff[m] = 0.0;
				for(int n1=0; n1<5; n1++){
					uf[m] += UI[n1][i][j][m+k-stencil_P]*eigen_l[n][n1];
					ff[m] += F_z[n1][i][j][m+k-stencil_P]*eigen_l[n][n1];
				}
				pp[m] = 0.5*(ff[m] + artificial_viscosity*uf[m]);  
				mm[m] = 0.5*(ff[m] - artificial_viscosity*uf[m]);  
			}
	
			// calculate the scalar numerical flux at x direction
			f_flux = (*CurrentInterpolation_P)(&pp[stencil_P], dz) 
				   + (*CurrentInterpolation_M)(&mm[stencil_P], dz);
		}
			// get Fp
			for(int n1=0; n1<5; n1++)
			_p[n][n1] = f_flux*eigen_r[n1][n];
	}
	// reconstruction the F_z flux terms
	for(int n=0; n<5; n++){
		F_z_wall[n][i][j][k] = 0.0;
		for(int n1=0; n1<5; n1++) 
			F_z_wall[n][i][j][k] += _p[n1][n];
	}
}
//-------------------------------------------------------------------------------------------------
//		dissipative flux with sixth order: no visocus effects
//-------------------------------------------------------------------------------------------------
void BaseFluid::PhysicalDissipation_0()
{

}
//-------------------------------------------------------------------------------------------------
//		dissipative flux with sixth order: with visocus effects
//-------------------------------------------------------------------------------------------------
void BaseFluid::PhysicalDissipation_1()
{
	double ducx[Xmax][Ymax][Zmax], dvcy[Xmax][Ymax][Zmax], dwcz[Xmax][Ymax][Zmax];
	double ducy[Xmax][Ymax][Zmax], dvcx[Xmax][Ymax][Zmax], dwcx[Xmax][Ymax][Zmax];
	double ducz[Xmax][Ymax][Zmax], dvcz[Xmax][Ymax][Zmax], dwcy[Xmax][Ymax][Zmax];
	//calculate cell center dirivatives in flow domian
	for(int i=4; i<=l+4; i++)
		for(int j=4; j<=d+4; j++)
			for(int k=4; k<=width+4; k++){
				ducx[i][j][k] = (45.0*(u[i+1][j][k] - u[i-1][j][k]) - 9.0*(u[i+2][j][k] - u[i-2][j][k]) 
							  + (u[i+3][j][k] - u[i-3][j][k]))/dx/60.0;
				dvcx[i][j][k] = (45.0*(v[i+1][j][k] - v[i-1][j][k]) - 9.0*(v[i+2][j][k] - v[i-2][j][k]) 
							  + (v[i+3][j][k] - v[i-3][j][k]))/dx/60.0;
				dwcx[i][j][k] = (45.0*(w[i+1][j][k] - w[i-1][j][k]) - 9.0*(w[i+2][j][k] - w[i-2][j][k]) 
							  + (w[i+3][j][k] - w[i-3][j][k]))/dx/60.0;
				ducy[i][j][k] = (45.0*(u[i][j+1][k] - u[i][j-1][k]) - 9.0*(u[i][j+2][k] - u[i][j-2][k]) 
							  + (u[i][j+3][k] - u[i][j-3][k]))/dy/60.0;
				dvcy[i][j][k] = (45.0*(v[i][j+1][k] - v[i][j-1][k]) - 9.0*(v[i][j+2][k] - v[i][j-2][k]) 
							  + (v[i][j+3][k] - v[i][j-3][k]))/dy/60.0;
				dwcy[i][j][k] = (45.0*(w[i][j+1][k] - w[i][j-1][k]) - 9.0*(w[i][j+2][k] - w[i][j-2][k]) 
							  + (w[i][j+3][k] - w[i][j-3][k]))/dy/60.0;
				ducz[i][j][k] = (45.0*(u[i][j][k+1] - u[i][j][k-1]) - 9.0*(u[i][j][k+2] - u[i][j][k-2]) 
							  + (u[i][j][k+3] - u[i][j][k-3]))/dz/60.0;
				dvcz[i][j][k] = (45.0*(v[i][j][k+1] - v[i][j][k-1]) - 9.0*(v[i][j][k+2] - v[i][j][k-2]) 
							  + (v[i][j][k+3] - v[i][j][k-3]))/dz/60.0;
				dwcz[i][j][k] = (45.0*(w[i][j][k+1] - w[i][j][k-1]) - 9.0*(w[i][j][k+2] - w[i][j][k-2]) 
							  + (w[i][j][k+3] - w[i][j][k-3]))/dz/60.0;
			}

	//boundary condition for derivatives
	//x direction
	for(int k=4; k<=width+4; k++) 
		for(int j=4; j<=d+4; j++) {
		//for left side
		for(int i=3; i>=0; i--){
			switch(xBl) {
			//the rigid wall conditions 	
			case 0:
				ducy[i][j][k] = ducy[7-i][j][k];
				ducz[i][j][k] = ducz[7-i][j][k];
				dvcy[i][j][k] = dvcy[7-i][j][k];
				dwcz[i][j][k] = dwcz[7-i][j][k];
				break;
			//the peroidic conditions	
			case 1:
				ducy[i][j][k] = ducy[l+i+1][j][k];
				ducz[i][j][k] = ducz[l+i+1][j][k];
				dvcy[i][j][k] = dvcy[l+i+1][j][k];
				dwcz[i][j][k] = dwcz[l+i+1][j][k];
				break;
			//the inflow conditions 	
			case 2:
				ducy[i][j][k] = 0.0;
				ducz[i][j][k] = 0.0;
				dvcy[i][j][k] = 0.0;
				dwcz[i][j][k] = 0.0;
				break;
			//the outflow conditions	
			case 3:
				ducy[i][j][k] = ducy[i+1][j][k];;
				ducz[i][j][k] = ducz[i+1][j][k];
				dvcy[i][j][k] = dvcy[i+1][j][k];;
				dwcz[i][j][k] = dwcz[i+1][j][k];
				break;
			//the no-slip conditions	
			case 4:
				ducy[i][j][k] = -ducy[7-i][j][k];
				ducz[i][j][k] = -ducz[7-i][j][k];
				dvcy[i][j][k] = -dvcy[7-i][j][k];
				dwcz[i][j][k] = -dwcz[7-i][j][k];
				break;
			}//end-switch
		}
		
		//for right side
		for(int i=l+5; i<=l+8; i++){
			switch(xBr) {
			//the rigid wall conditions 	
			case 0:
				ducy[i][j][k] = ducy[2*l+9-i][j][k];
				ducz[i][j][k] = ducz[2*l+9-i][j][k];
				dvcy[i][j][k] = dvcy[2*l+9-i][j][k];
				dwcz[i][j][k] = dwcz[2*l+9-i][j][k];
				break;
			//the peroidic conditions	
			case 1:
				ducy[i][j][k] = ducy[i-l-1][j][k];
				ducz[i][j][k] = ducz[i-l-1][j][k];
				dvcy[i][j][k] = dvcy[i-l-1][j][k];
				dwcz[i][j][k] = dwcz[i-l-1][j][k];
				break;
			//the inflow conditions 	
			case 2:
				ducy[i][j][k] = 0.0;
				ducz[i][j][k] = 0.0;
				dvcy[i][j][k] = 0.0;
				dwcz[i][j][k] = 0.0;
				break;
			//the outflow conditions	
			case 3:
				ducy[i][j][k] = ducy[i-1][j][k];
				ducz[i][j][k] = ducz[i-1][j][k];
				dvcy[i][j][k] = dvcy[i-1][j][k];
				dwcz[i][j][k] = dwcz[i-1][j][k];
				break;
			//the no-slip conditions 	
			case 4:
				ducy[i][j][k] = -ducy[2*l+9-i][j][k];
				ducz[i][j][k] = -ducz[2*l+9-i][j][k];
				dvcy[i][j][k] = -dvcy[2*l+9-i][j][k];
				dwcz[i][j][k] = -dwcz[2*l+9-i][j][k];
				break;
			}//end-switch
		}
	}

	//y direction
	for(int k=4; k<=width+4; k++)
		for(int i=0; i<=l+8; i++){
		//the bottom side
		for(int j=3; j>=0; j--) {
			switch(yBd) {
			//the rigid wall conditions 	
			case 0:
				dvcx[i][j][k] = dvcx[i][7-j][k];
				dvcz[i][j][k] = dvcz[i][7-j][k];
				ducx[i][j][k] = ducx[i][7-j][k];
				dwcz[i][j][k] = dwcz[i][7-j][k];
				break;
			//the peroidic conditions	
			case 1:
				dvcx[i][j][k] = dvcx[i][j+d+1][k];
				dvcz[i][j][k] = dvcz[i][j+d+1][k];
				ducx[i][j][k] = ducx[i][j+d+1][k];
				dwcz[i][j][k] = dwcz[i][j+d+1][k];
				break;
			//the inflow conditions 	
			case 2:
				dvcx[i][j][k] = 0.0;
				dvcz[i][j][k] = 0.0;
				ducx[i][j][k] = 0.0;
				dwcz[i][j][k] = 0.0;
				break;
				//the outflow conditions	
			case 3:
				dvcx[i][j][k] = dvcx[i][j+1][k];
				dvcz[i][j][k] = dvcz[i][j+1][k];
				ducx[i][j][k] = ducx[i][j+1][k];
				dwcz[i][j][k] = dwcz[i][j+1][k];
				break;
			//the no-slip conditions 	
			case 4:
				dvcx[i][j][k] = -dvcx[i][7-j][k];
				dvcz[i][j][k] = -dvcz[i][7-j][k];
				ducx[i][j][k] = -ducx[i][7-j][k];
				dwcz[i][j][k] = -dwcz[i][7-j][k];
			break;
			}//end-switch
		}
				
			//the top side
			for(int j=d+5; j<=d+8; j++) {
				switch(yBu) {
				//the rigid wall conditions 	
				case 0:
					dvcx[i][j][k] = dvcx[i][2*d+9-j][k];
					dvcz[i][j][k] = dvcz[i][2*d+9-j][k];
					ducx[i][j][k] = ducx[i][2*d+9-j][k];
					dwcz[i][j][k] = dwcz[i][2*d+9-j][k];
					break;
					//the peroidic conditions	
				case 1:
					dvcx[i][j][k] = dvcx[i][j-d-1][k];
					dvcz[i][j][k] = dvcz[i][j-d-1][k];
					ducx[i][j][k] = ducx[i][j-d-1][k];
					dwcz[i][j][k] = dwcz[i][j-d-1][k];
					break;
				//the inflow conditions 	
				case 2:
					dvcx[i][j][k] = 0.0;
					dvcz[i][j][k] = 0.0;
					ducx[i][j][k] = 0.0;
					dwcz[i][j][k] = 0.0;
					break;
				//the outflow conditions	
				case 3:
					dvcx[i][j][k] = dvcx[i][j-1][k];
					dvcz[i][j][k] = dvcz[i][j-1][k];
					ducx[i][j][k] = ducx[i][j-1][k];
					dwcz[i][j][k] = dwcz[i][j-1][k];
					break;
				//the no-slip conditions 	
				case 4:
					dvcx[i][j][k] = -dvcx[i][2*d+9-j][k];
					dvcz[i][j][k] = -dvcz[i][2*d+9-j][k];
					ducx[i][j][k] = -ducx[i][2*d+9-j][k];
					dwcz[i][j][k] = -dwcz[i][2*d+9-j][k];
					break;
				}//end-switch
			}
		}
	
	//z direction
	for(int j=0; j<=d+8; j++)
		for(int i=0; i<=l+8; i++){
			//the near side
			for(int k=3; k>=0; k--) {
				switch(zBn) {
				//the rigid wall conditions 	
				case 0:
					dwcx[i][j][k] = dwcx[i][j][7-k];
					dwcy[i][j][k] = dwcy[i][j][7-k];
					ducx[i][j][k] = ducx[i][j][7-k];
					dvcy[i][j][k] = dvcy[i][j][7-k];
					break;
				//the peroidic conditions	
				case 1:
					dwcx[i][j][k] = dwcx[i][j][k+width+1];
					dwcy[i][j][k] = dwcy[i][j][k+width+1];
					ducx[i][j][k] = ducx[i][j][k+width+1];
					dvcy[i][j][k] = dvcy[i][j][k+width+1];
					break;
				//the inflow conditions 	
				case 2:
					dwcx[i][j][k] = 0.0;
					dwcy[i][j][k] = 0.0;
					ducx[i][j][k] = 0.0;
					dvcy[i][j][k] = 0.0;
					break;
					//the outflow conditions	
				case 3:
					dwcx[i][j][k] = dwcx[i][j][k+1];
					dwcy[i][j][k] = dwcy[i][j][k+1];
					ducx[i][j][k] = ducx[i][j][k+1];
					dvcy[i][j][k] = dvcy[i][j][k+1];
					break;
				//the no-slip conditions 	
				case 4:
					dwcx[i][j][k] = -dwcx[i][j][7-k];
					dwcy[i][j][k] = -dwcy[i][j][7-k];
					ducx[i][j][k] = -ducx[i][j][7-k];
					dvcy[i][j][k] = -dvcy[i][j][7-k];
					break;
					}//end-switch
			}
				
			//the far side
			for(int k=width+5; k<=width+8; k++) {
				switch(zBf) {
				//the rigid wall conditions 	
				case 0:
					dwcx[i][j][k] = dwcx[i][j][2*width+9-k];
					dwcy[i][j][k] = dwcy[i][j][2*width+9-k];
					ducx[i][j][k] = ducx[i][j][2*width+9-k];
					dvcy[i][j][k] = dvcy[i][j][2*width+9-k];
					break;
					//the peroidic conditions	
				case 1:
					dwcx[i][j][k] = dwcx[i][j][k-width-1];
					dwcy[i][j][k] = dwcy[i][j][k-width-1];
					ducx[i][j][k] = ducx[i][j][k-width-1];
					dvcy[i][j][k] = dvcy[i][j][k-width-1];
					break;
				//the inflow conditions 	
				case 2:
					dwcx[i][j][k] = 0.0;
					dwcy[i][j][k] = 0.0;
					ducx[i][j][k] = 0.0;
					dvcy[i][j][k] = 0.0;
					break;
				//the outflow conditions	
				case 3:
					dwcx[i][j][k] = dwcx[i][j][k-1];
					dwcy[i][j][k] = dwcy[i][j][k-1];
					ducx[i][j][k] = ducx[i][j][k-1];
					dvcy[i][j][k] = dvcy[i][j][k-1];
					break;
				//the no-slip conditions 	
				case 4:
					dwcx[i][j][k] = -dwcx[i][j][2*width+9-k];
					dwcy[i][j][k] = -dwcy[i][j][2*width+9-k];
					ducx[i][j][k] = -ducx[i][j][2*width+9-k];
					dvcy[i][j][k] = -dvcy[i][j][2*width+9-k];
					break;
				}//end-switch
			}
		}			

	//dissiaptive fluxes in x direction
	for(int i=3; i<=l+4; i++)
		for(int j=4; j<=d+4; j++)
			for(int k=4; k<=width+4; k++){
				double f_x = (2.0*mu1 + mu2)*(2250.0*(u[i+1][j][k] - u[i][j][k]) - 125.0*(u[i+2][j][k] - u[i-1][j][k]) 
						   + 9.0*(u[i+3][j][k] - u[i-2][j][k]))/dx/1920.0;
				double f_y = mu1*(2250.0*(v[i+1][j][k] - v[i][j][k]) - 125.0*(v[i+2][j][k] - v[i-1][j][k]) 
						   + 9.0*(v[i+3][j][k] - v[i-2][j][k]))/dx/1920.0;
				double f_z = mu1*(2250.0*(w[i+1][j][k] - w[i][j][k]) - 125.0*(w[i+2][j][k] - w[i-1][j][k]) 
						   + 9.0*(w[i+3][j][k] - w[i-2][j][k]))/dx/1920.0;
				f_x += mu2*(150.0*(dvcy[i+1][j][k] + dvcy[i][j][k]) - 25.0*(dvcy[i+2][j][k] + dvcy[i-1][j][k])
					+ 3.0*(dvcy[i+3][j][k] + dvcy[i-2][j][k])  + 150.0*(dwcz[i+1][j][k] + dwcz[i][j][k]) 
					- 25.0*(dwcz[i+2][j][k] + dwcz[i-1][j][k]) + 3.0*(dwcz[i+3][j][k] + dwcz[i-2][j][k]))/256.0;
				f_y += mu1*(150.0*(ducy[i+1][j][k] + ducy[i][j][k]) - 25.0*(ducy[i+2][j][k] + ducy[i-1][j][k])
							  + 3.0*(ducy[i+3][j][k] + ducy[i-2][j][k]))/256.0;
				f_z += mu1*(150.0*(ducz[i+1][j][k] + ducz[i][j][k]) - 25.0*(ducz[i+2][j][k] + ducz[i-1][j][k])
							  + 3.0*(ducz[i+3][j][k] + ducz[i-2][j][k]))/256.0;
				double u_hlf = (150.0*(u[i+1][j][k] + u[i][j][k]) - 25.0*(u[i+2][j][k] + u[i-1][j][k])
							  + 3.0*(u[i+3][j][k] + u[i-2][j][k]))/256.0;
				double v_hlf = (150.0*(v[i+1][j][k] + v[i][j][k]) - 25.0*(v[i+2][j][k] + v[i-1][j][k])
							  + 3.0*(v[i+3][j][k] + v[i-2][j][k]))/256.0;
				double w_hlf = (150.0*(w[i+1][j][k] + w[i][j][k]) - 25.0*(w[i+2][j][k] + w[i-1][j][k])
							  + 3.0*(w[i+3][j][k] + w[i-2][j][k]))/256.0;
				F_x_wall[1][i][j][k] -= f_x;
				F_x_wall[2][i][j][k] -= f_y;
				F_x_wall[3][i][j][k] -= f_z;
				F_x_wall[4][i][j][k] -= f_x*u_hlf + f_y*v_hlf + f_z*w_hlf;
				F_x_wall[4][i][j][k] -= rho_alpha*(2250.0*(e[i+1][j][k] - e[i][j][k]) - 125.0*(e[i+2][j][k] - e[i-1][j][k]) 
									  + 9.0*(e[i+3][j][k] - e[i-2][j][k]))/dx/1920.0;;
			}

	if(d != 0) {
		//dissiaptive fluxes in y direction
		for(int i=4; i<=l+4; i++)
			for(int k=4; k<=width+4; k++)
				for(int j=3; j<=d+4; j++){
					double f_x = mu1*(2250.0*(u[i][j+1][k] - u[i][j][k]) - 125.0*(u[i][j+2][k] - u[i][j-1][k]) 
							   + 9.0*(u[i][j+3][k] - u[i][j-2][k]))/dy/1920.0;
					double f_y = (2.0*mu1 + mu2)*(2250.0*(v[i][j+1][k] - v[i][j][k]) - 125.0*(v[i][j+2][k] - v[i][j-1][k]) 
							   + 9.0*(v[i][j+3][k] - v[i][j-2][k]))/dy/1920.0;
					double f_z = mu1*(2250.0*(w[i][j+1][k] - w[i][j][k]) - 125.0*(w[i][j+2][k] - w[i][j-1][k]) 
							   + 9.0*(w[i][j+3][k] - w[i][j-2][k]))/dy/1920.0;
					f_x += mu1*(150.0*(dvcx[i][j+1][k] + dvcx[i][j][k]) - 25.0*(dvcx[i][j+2][k] + dvcx[i][j-1][k])
								  + 3.0*(dvcx[i][j+3][k] + dvcx[i][j-2][k]))/256.0;
					f_y += mu2*(150.0*(ducx[i][j+1][k] + ducx[i][j][k]) - 25.0*(ducx[i][j+2][k] + ducx[i][j-1][k])
						+ 3.0*(ducx[i][j+3][k] + ducx[i][j-2][k])  + 150.0*(dwcz[i][j+1][k] + dwcz[i][j][k]) 
						- 25.0*(dwcz[i][j+2][k] + dwcz[i][j-1][k]) + 3.0*(dwcz[i][j+3][k] + dwcz[i][j-2][k]))/256.0/3.0;
					f_z += mu1*(150.0*(dvcz[i][j+1][k] + dvcz[i][j][k]) - 25.0*(dvcz[i][j+2][k] + dvcz[i][j-1][k])
								  + 3.0*(dvcz[i][j+3][k] + dvcz[i][j-2][k]))/256.0;
					double u_hlf = (150.0*(u[i][j+1][k] + u[i][j][k]) - 25.0*(u[i][j+2][k] + u[i][j-1][k])
								  + 3.0*(u[i][j+3][k] + u[i][j-2][k]))/256.0;
					double v_hlf = (150.0*(v[i][j+1][k] + v[i][j][k]) - 25.0*(v[i][j+2][k] + v[i][j-1][k])
								  + 3.0*(v[i][j+3][k] + v[i][j-2][k]))/256.0;
					double w_hlf = (150.0*(w[i][j+1][k] + w[i][j][k]) - 25.0*(w[i][j+2][k] + w[i][j-1][k])
								  + 3.0*(w[i][j+3][k] + w[i][j-2][k]))/256.0;
					F_y_wall[1][i][j][k] -= f_x;
					F_y_wall[2][i][j][k] -= f_y;
					F_y_wall[3][i][j][k] -= f_z;
					F_y_wall[4][i][j][k] -= f_x*u_hlf + f_y*v_hlf + f_z*w_hlf;
					F_y_wall[4][i][j][k] -= rho_alpha*(2250.0*(e[i][j+1][k] - e[i][j][k]) - 125.0*(e[i][j+2][k] - e[i][j-1][k]) 
										  + 9.0*(e[i][j+3][k] - e[i][j-2][k]))/dy/1920.0;
				}
	}

	if(width != 0) {
		//dissiaptive fluxes in z direction
		for(int i=4; i<=l+4; i++)
			for(int j=4; j<=d+4; j++)
				//get the G-flux terms at cell wall
				for(int k=3; k<=width+4; k++){
					double f_x = mu1*(2250.0*(u[i][j][k+1] - u[i][j][k]) - 125.0*(u[i][j][k+2] - u[i][j][k-1]) 
							   + 9.0*(u[i][j][k+3] - u[i][j][k-2]))/dz/1920.0;
					double f_y = mu1*(2250.0*(v[i][j][k+1] - v[i][j][k]) - 125.0*(v[i][j][k+2] - v[i][j][k-1]) 
							   + 9.0*(v[i][j][k+3] - v[i][j][k-2]))/dz/1920.0;
					double f_z = (2.0*mu1 + mu2)*(2250.0*(w[i][j][k+1] - w[i][j][k]) - 125.0*(w[i][j][k+2] - w[i][j][k-1]) 
							   + 9.0*(w[i][j][k+3] - w[i][j][k-2]))/dz/1920.0;
					f_x += mu1*(150.0*(dwcx[i][j][k+1] + dwcx[i][j][k]) - 25.0*(dwcx[i][j][k+2] + dwcx[i][j][k-1])
								  + 3.0*(dwcx[i][j][k+3] + dwcx[i][j][k-2]))/256.0;
					f_y += mu1*(150.0*(dwcy[i][j][k+1] + dwcy[i][j][k]) - 25.0*(dwcy[i][j][k+2] + dwcy[i][j][k-1])
								  + 3.0*(dwcy[i][j][k+3] + dwcy[i][j][k-2]))/256.0;
					f_z += mu2*(150.0*(ducx[i][j+1][k] + ducx[i][j][k]) - 25.0*(ducx[i][j+2][k] + ducx[i][j-1][k])
						+ 3.0*(ducx[i][j+3][k] + ducx[i][j-2][k])  + 150.0*(dvcy[i][j+1][k] + dvcy[i][j][k]) 
						- 25.0*(dvcy[i][j+2][k] + dvcy[i][j-1][k]) + 3.0*(dvcy[i][j+3][k] + dvcy[i][j-2][k]))/256.0/3.0;
					double u_hlf = (150.0*(u[i][j][k+1] + u[i][j][k]) - 25.0*(u[i][j][k+2] + u[i][j][k-1])
								  + 3.0*(u[i][j][k+3] + u[i][j][k-2]))/256.0;
					double v_hlf = (150.0*(v[i][j][k+1] + v[i][j][k]) - 25.0*(v[i][j][k+2] + v[i][j][k-1])
								  + 3.0*(v[i][j][k+3] + v[i][j][k-2]))/256.0;
					double w_hlf = (150.0*(w[i][j][k+1] + w[i][j][k]) - 25.0*(w[i][j][k+2] + w[i][j][k-1])
								  + 3.0*(w[i][j][k+3] + w[i][j][k-2]))/256.0;				
					F_z_wall[1][i][j][k] -= f_x;
					F_z_wall[2][i][j][k] -= f_y;
					F_z_wall[3][i][j][k] -= f_z;
					F_z_wall[4][i][j][k] -= f_x*u_hlf + f_y*v_hlf + f_z*w_hlf;
					F_z_wall[4][i][j][k] -= rho_alpha*(2250.0*(e[i][j][k+1] - e[i][j][k]) - 125.0*(e[i][j][k+2] - e[i][j][k-1]) 
										  + 9.0*(e[i][j][k+3] - e[i][j][k-2]))/dz/1920.0;
				}
	}
}
//-------------------------------------------------------------------------------------------------
//		Positivity preserving limiter not active
//-------------------------------------------------------------------------------------------------
void BaseFluid::PositivityPreserving_0(double UI[Emax][Xmax][Ymax][Zmax])
{

}
//-------------------------------------------------------------------------------------------------
//		Positivity preserving limiter active
//-------------------------------------------------------------------------------------------------
void BaseFluid::PositivityPreserving_1(double UI[Emax][Xmax][Ymax][Zmax])
{
	double epsilon = 1.0e-13;
	//proceed at x directiom and get F-flux terms at node wall
	for(int j=4; j<=d+4; j++)
		for(int k=4; k<=width+4; k++)
			//get the F-flux terms at cell wall
			for(int i=3; i<=l+4; i++){
				//a partition of conservative variables
				double rho_min, p_min, theta, theta_u, theta_p;
				double UU[5], UP[5], F_LF[5], FF_LF[5], FF[5];
				for(int n=0; n<5; n++) {
					UU[n] = UI[n][i][j][k];
					UP[n] = UI[n][i+1][j][k];
					F_LF[n] = 0.5*(F_x[n][i][j][k] + F_x[n][i+1][j][k] + lambda_x0*(UI[n][i][j][k] - UI[n][i+1][j][k]));
					FF_LF[n] = 2.0*lambda_x*F_LF[n];
					FF[n] = 2.0*lambda_x*F_x_wall[n][i][j][k];
				}
				//correct for positive density
				theta_u = 1.0; theta_p = 1.0;
				rho_min = AMIN1(UU[0], epsilon);
				if(UU[0] - FF[0] < rho_min) theta_u = (UU[0] - FF_LF[0] - rho_min)/(FF[0] - FF_LF[0]);
				rho_min =  AMIN1(UP[0], epsilon);
				if(UP[0] + FF[0] < rho_min)  theta_p = (UP[0] + FF_LF[0] - rho_min)/(FF_LF[0] - FF[0]);
				theta = AMIN1(theta_u, theta_p);
				for(int n=0; n<5; n++) {
					FF[n] = (1.0 - theta)*FF_LF[n] + theta*FF[n]; 
					F_x_wall[n][i][j][k] = (1.0 - theta)*F_LF[n] + theta*F_x_wall[n][i][j][k];
				}
				//correct for positive pressure
				theta_u = 1.0; theta_p = 1.0;
				double p_q = _Gamma*(UU[4] - FF[4] - 0.5*((UU[1] - FF[1])*(UU[1] - FF[1])
							+ (UU[2] - FF[2])*(UU[2] - FF[2]) 
							+ (UU[3] - FF[3])*(UU[3] - FF[3]))/(UU[0] - FF[0]));
				if(p_q < epsilon) {
					double p_u = _Gamma*(UU[4] - FF_LF[4] - 0.5*((UU[1] - FF_LF[1])*(UU[1] - FF_LF[1])
								+ (UU[2] - FF_LF[2])*(UU[2] - FF_LF[2]) 
								+ (UU[3] - FF_LF[3])*(UU[3] - FF_LF[3]))/(UU[0] - FF_LF[0]));
					p_min =  AMIN1(p_u, epsilon);
					theta_u = (p_u - p_min)/(p_u - p_q);
				}
				p_q = _Gamma*(UP[4] + FF[4] - 0.5*((UP[1] + FF[1])*(UP[1] + FF[1])
							+ (UP[2] + FF[2])*(UP[2] + FF[2]) 
							+ (UP[3] + FF[3])*(UP[3] + FF[3]))/(UP[0] + FF[0]));
				if(p_q < epsilon) {
					double p_p = _Gamma*(UP[4] + FF_LF[4] - 0.5*((UP[1] + FF_LF[1])*(UP[1] + FF_LF[1])
							+ (UP[2] + FF_LF[2])*(UP[2] + FF_LF[2]) 
							+ (UP[3] + FF_LF[3])*(UP[3] + FF_LF[3]))/(UP[0] + FF_LF[0]));
					p_min =  AMIN1(p_p, epsilon);
					theta_p = (p_p - p_min)/(p_p - p_q);
				}
				theta = AMIN1(theta_u, theta_p);
				for(int n=0; n<5; n++) F_x_wall[n][i][j][k] = (1.0 - theta)*F_LF[n] + theta*F_x_wall[n][i][j][k];
	}
		
	//proceed at y directiom and get G-flux terms at node wall
	for(int i=4; i<=l+4; i++)
		for(int k=4; k<=width+4; k++)
			//get the G-flux terms at cell wall
			for(int j=3; j<=d+4; j++){
				//a partition of conservative variables
				double rho_min, p_min, theta, theta_u, theta_p;
				double UU[5], UP[5], F_LF[5], FF_LF[5], FF[5];
				for(int n=0; n<5; n++) {
					UU[n] = UI[n][i][j][k];
					UP[n] = UI[n][i][j+1][k];
					F_LF[n] = 0.5*(F_y[n][i][j][k] + F_y[n][i][j+1][k] + lambda_y0*(UI[n][i][j][k] - UI[n][i][j+1][k]));
					FF_LF[n] = 2.0*lambda_y*F_LF[n];
					FF[n] = 2.0*lambda_y*F_y_wall[n][i][j][k];
				}
				//correct for positive density
				theta_u = 1.0; theta_p = 1.0;
				rho_min = AMIN1(UU[0], epsilon);
				if(UU[0] - FF[0] < rho_min) theta_u = (UU[0] - FF_LF[0] - rho_min)/(FF[0] - FF_LF[0]);
				rho_min = AMIN1(UP[0], epsilon);
				if(UP[0] + FF[0] < rho_min) theta_p = (UP[0] + FF_LF[0] - rho_min)/(FF_LF[0] - FF[0]);
				theta = AMIN1(theta_u, theta_p);
				for(int n=0; n<5; n++) {
					FF[n] = (1.0 - theta)*FF_LF[n] + theta*FF[n]; 
					F_y_wall[n][i][j][k] = (1.0 - theta)*F_LF[n] + theta*F_y_wall[n][i][j][k];
				}
				//correct for positive pressure
				theta_u = 1.0; theta_p = 1.0;
				double p_q = _Gamma*(UU[4] - FF[4] - 0.5*((UU[1] - FF[1])*(UU[1] - FF[1])
							+ (UU[2] - FF[2])*(UU[2] - FF[2]) 
							+ (UU[3] - FF[3])*(UU[3] - FF[3]))/(UU[0] - FF[0]));
				if(p_q < epsilon) {
					double p_u = _Gamma*(UU[4] - FF_LF[4] - 0.5*((UU[1] - FF_LF[1])*(UU[1] - FF_LF[1])
								+ (UU[2] - FF_LF[2])*(UU[2] - FF_LF[2]) 
								+ (UU[3] - FF_LF[3])*(UU[3] - FF_LF[3]))/(UU[0] - FF_LF[0]));
					p_min = AMIN1(p_u, epsilon);
					theta_u = (p_u - p_min)/(p_u - p_q);
				}
				p_q = _Gamma*(UP[4] + FF[4] - 0.5*((UP[1] + FF[1])*(UP[1] + FF[1])
							+ (UP[2] + FF[2])*(UP[2] + FF[2]) 
							+ (UP[3] + FF[3])*(UP[3] + FF[3]))/(UP[0] + FF[0]));
				if(p_q < epsilon) {
					double p_p = _Gamma*(UP[4] + FF_LF[4] - 0.5*((UP[1] + FF_LF[1])*(UP[1] + FF_LF[1])
							+ (UP[2] + FF_LF[2])*(UP[2] + FF_LF[2]) 
							+ (UP[3] + FF_LF[3])*(UP[3] + FF_LF[3]))/(UP[0] + FF_LF[0]));
					p_min = AMIN1(p_p, epsilon);
					theta_p = (p_p - p_min)/(p_p - p_q);
				}
				theta = AMIN1(theta_u, theta_p);
				for(int n=0; n<5; n++) F_y_wall[n][i][j][k] = (1.0 - theta)*F_LF[n] + theta*F_y_wall[n][i][j][k];
	}

	//proceed at z directiom and get G-flux terms at node wall
//	for(int i=4; i<=l+4; i++)
//		for(int j=4; j<=d+4; j++)
			//get the G-flux terms at cell wall
//			for(int k=3; k<=width+4; k++){
//	}	
}
//-------------------------------------------------------------------------------------------------
//		get Roe average
//-------------------------------------------------------------------------------------------------
void BaseFluid::RoeAverage_x(int i, int j, int k, double eigen_l[Emax][Emax], 
							 double eigen_r[Emax][Emax], double eigen_value[Emax], double &ave_rho1) 
{
	//preparing some interval value
	double D	=	sqrt(rho[i+1][j][k]/rho[i][j][k]);
	double D1	=	1.0 / (D + 1.0);
	double _u	=	(u[i][j][k] + D*u[i+1][j][k])*D1;
	double _v	=	(v[i][j][k] + D*v[i+1][j][k])*D1;
	double _w	=	(w[i][j][k] + D*w[i+1][j][k])*D1;
	double _H	=	(H[i][j][k] + D*H[i+1][j][k])*D1;
	double _rho = sqrt(rho[i+1][j][k]*rho[i][j][k]);
	double _rho1 = 1.0 / _rho;
	double q2 = _u*_u + _v*_v + _w*_w;
	double c2 = _Gamma*(_H - 0.5*q2);
	double _c = sqrt(c2);
	double _c1_rho = 0.5*_rho / _c;
	double c21_Gamma = _Gamma / c2;
	double _c1_rho1_Gamma = _Gamma*_rho1 / _c;

	// left eigen vectors 
	eigen_l[0][0] = 1.0 - 0.5*c21_Gamma*q2;
	eigen_l[0][1] = c21_Gamma*_u;
	eigen_l[0][2] = c21_Gamma*_v;
	eigen_l[0][3] = c21_Gamma*_w;
	eigen_l[0][4] = -c21_Gamma;
	
	eigen_l[1][0] = -_w*_rho1;
	eigen_l[1][1] = 0.0;
	eigen_l[1][2] = 0.0;
	eigen_l[1][3] = _rho1;
	eigen_l[1][4] = 0.0;

	eigen_l[2][0] = _v*_rho1;
	eigen_l[2][1] = 0.0;
	eigen_l[2][2] = -_rho1;
	eigen_l[2][3] = 0.0;
	eigen_l[2][4] = 0.0;

	eigen_l[3][0] = 0.5*_c1_rho1_Gamma*q2 - _u*_rho1;
	eigen_l[3][1] = -_c1_rho1_Gamma*_u + _rho1;
	eigen_l[3][2] = -_c1_rho1_Gamma*_v;
	eigen_l[3][3] = -_c1_rho1_Gamma*_w;
	eigen_l[3][4] = _c1_rho1_Gamma;

	eigen_l[4][0] = 0.5*_c1_rho1_Gamma*q2 + _u*_rho1;
	eigen_l[4][1] = -_c1_rho1_Gamma*_u - _rho1;
	eigen_l[4][2] = -_c1_rho1_Gamma*_v;
	eigen_l[4][3] = -_c1_rho1_Gamma*_w;
	eigen_l[4][4] = _c1_rho1_Gamma;

	//right eigen vectors
	eigen_r[0][0] = 1.0;
	eigen_r[0][1] = 0.0;
	eigen_r[0][2] = 0.0;
	eigen_r[0][3] = _c1_rho;
	eigen_r[0][4] = _c1_rho;
	
	eigen_r[1][0] = _u;
	eigen_r[1][1] = 0.0;
	eigen_r[1][2] = 0.0;
	eigen_r[1][3] = _c1_rho*(_u + _c);
	eigen_r[1][4] = _c1_rho*(_u - _c);

	eigen_r[2][0] = _v;
	eigen_r[2][1] = 0.0;
	eigen_r[2][2] = -_rho;
	eigen_r[2][3] = _c1_rho*_v;
	eigen_r[2][4] = _c1_rho*_v;

	eigen_r[3][0] = _w;
	eigen_r[3][1] = _rho;
	eigen_r[3][2] = 0.0;
	eigen_r[3][3] = _c1_rho*_w;
	eigen_r[3][4] = _c1_rho*_w;

	eigen_r[4][0] = 0.5*q2;
	eigen_r[4][1] = _rho*_w;
	eigen_r[4][2] = -_rho*_v;
	eigen_r[4][3] = _c1_rho*(_H + _u*_c);
	eigen_r[4][4] = _c1_rho*(_H - _u*_c);

	eigen_value[0] = fabs(_u);
	eigen_value[1] = eigen_value[0];
	eigen_value[2] = eigen_value[0];
//	double du = fabs(u[i][j][k] - u[i + 1][j][k]);
//	_c = AMIN1(0.0*du, _c);
	eigen_value[3] = fabs(_u + _c);
	eigen_value[4] = fabs(_u - _c);
	ave_rho1 = _rho1;
}
void BaseFluid::RoeAverage_y(int i, int j, int k, double eigen_l[Emax][Emax], 
							 double eigen_r[Emax][Emax], double eigen_value[Emax], double &ave_rho1) 
{
	//preparing some interval value
	double D	=	sqrt(rho[i][j+1][k]/rho[i][j][k]);
	double D1	=	1.0 / (D + 1.0);
	double _u	=	(u[i][j][k] + D*u[i][j+1][k])*D1;
	double _v	=	(v[i][j][k] + D*v[i][j+1][k])*D1;
	double _w	=	(w[i][j][k] + D*w[i][j+1][k])*D1;
	double _H	=	(H[i][j][k] + D*H[i][j+1][k])*D1;
	double _rho = sqrt(rho[i][j+1][k]*rho[i][j][k]);
	double _rho1 = 1.0 / _rho;
	double q2 = _u*_u + _v*_v + _w*_w;
	double c2 = _Gamma*(_H - 0.5*q2);
	double _c = sqrt(c2);
	double _c1_rho = 0.5*_rho / _c;
	double c21_Gamma = _Gamma / c2;
	double _c1_rho1_Gamma = _Gamma*_rho1 / _c;

	// left eigen vectors 
	eigen_l[0][0] = _w*_rho1;
	eigen_l[0][1] = 0.0;
	eigen_l[0][2] = 0.0;
	eigen_l[0][3] = - _rho1;
	eigen_l[0][4] = 0.0;
	
	eigen_l[1][0] = 1.0 - 0.5*c21_Gamma*q2;
	eigen_l[1][1] = c21_Gamma*_u;
	eigen_l[1][2] = c21_Gamma*_v;
	eigen_l[1][3] = c21_Gamma*_w;
	eigen_l[1][4] = - c21_Gamma;

	eigen_l[2][0] = - _u*_rho1;
	eigen_l[2][1] = _rho1;
	eigen_l[2][2] = 0.0;
	eigen_l[2][3] = 0.0;
	eigen_l[2][4] = 0.0;

	eigen_l[3][0] = 0.5*_c1_rho1_Gamma*q2 - _v*_rho1;
	eigen_l[3][1] = -_c1_rho1_Gamma*_u;
	eigen_l[3][2] = -_c1_rho1_Gamma*_v + _rho1;
	eigen_l[3][3] = -_c1_rho1_Gamma*_w;
	eigen_l[3][4] = _c1_rho1_Gamma;

	eigen_l[4][0] = 0.5*_c1_rho1_Gamma*q2 + _v*_rho1;
	eigen_l[4][1] = -_c1_rho1_Gamma*_u;
	eigen_l[4][2] = -_c1_rho1_Gamma*_v - _rho1;
	eigen_l[4][3] = -_c1_rho1_Gamma*_w;
	eigen_l[4][4] = _c1_rho1_Gamma;

	//right eigen vectors
	eigen_r[0][0] = 0.0;
	eigen_r[0][1] = 1.0;
	eigen_r[0][2] = 0.0;
	eigen_r[0][3] = _c1_rho;
	eigen_r[0][4] = _c1_rho;
	
	eigen_r[1][0] = 0.0;
	eigen_r[1][1] = _u;
	eigen_r[1][2] = _rho;
	eigen_r[1][3] = _c1_rho*_u;
	eigen_r[1][4] = _c1_rho*_u;

	eigen_r[2][0] = 0.0;
	eigen_r[2][1] = _v;
	eigen_r[2][2] = 0.0;
	eigen_r[2][3] = _c1_rho*(_v + _c) ;
	eigen_r[2][4] = _c1_rho*(_v - _c);

	eigen_r[3][0] = - _rho;
	eigen_r[3][1] = _w;
	eigen_r[3][2] = 0.0;
	eigen_r[3][3] = _c1_rho*_w;
	eigen_r[3][4] = _c1_rho*_w;

	eigen_r[4][0] = - _rho*_w;
	eigen_r[4][1] = 0.5*q2;
	eigen_r[4][2] = _rho*_u;
	eigen_r[4][3] = _c1_rho*(_H + _v*_c);
	eigen_r[4][4] = _c1_rho*(_H - _v*_c);

	eigen_value[0] = fabs(_v);
	eigen_value[1] = eigen_value[0];
	eigen_value[2] = eigen_value[0];
	eigen_value[3] = fabs(_v + _c);
	eigen_value[4] = fabs(_v - _c);
	ave_rho1 = _rho1;
}
void BaseFluid::RoeAverage_z(int i, int j, int k, double eigen_l[Emax][Emax], 
							 double eigen_r[Emax][Emax], double eigen_value[Emax], double &ave_rho1)
{
	//preparing some interval value
	double D	=	sqrt(rho[i][j][k+1]/rho[i][j][k]);
	double D1	=	1.0 / (D + 1.0);
	double _u	=	(u[i][j][k] + D*u[i][j][k+1])*D1;
	double _v	=	(v[i][j][k] + D*v[i][j][k+1])*D1;
	double _w	=	(w[i][j][k] + D*w[i][j][k+1])*D1;
	double _H	=	(H[i][j][k] + D*H[i][j][k+1])*D1;
	double _rho = sqrt(rho[i][j][k+1]*rho[i][j][k]);
	double _rho1 = 1.0 / _rho;
	double q2 = _u*_u + _v*_v + _w*_w;
	double c2 = _Gamma*(_H - 0.5*q2);
	double _c = sqrt(c2);
	double _c1_rho = 0.5*_rho / _c;
	double c21_Gamma = _Gamma / c2;
	double _c1_rho1_Gamma = _Gamma*_rho1 / _c;

	// left eigen vectors 
	eigen_l[0][0] = - _v*_rho1;
	eigen_l[0][1] = 0.0;
	eigen_l[0][2] = _rho1;
	eigen_l[0][3] = 0.0;
	eigen_l[0][4] = 0.0;
	
	eigen_l[1][0] = _u*_rho1;
	eigen_l[1][1] = - _rho1;
	eigen_l[1][2] = 0.0;
	eigen_l[1][3] = 0.0;
	eigen_l[1][4] = 0.0;

	eigen_l[2][0] = 1.0 - 0.5*c21_Gamma*q2; 
	eigen_l[2][1] = c21_Gamma*_u; 
	eigen_l[2][2] = c21_Gamma*_v; 
	eigen_l[2][3] = c21_Gamma*_w; 
	eigen_l[2][4] = - c21_Gamma;

	eigen_l[3][0] = 0.5*_c1_rho1_Gamma*q2 - _w*_rho1;
	eigen_l[3][1] = -_c1_rho1_Gamma*_u;
	eigen_l[3][2] = -_c1_rho1_Gamma*_v;
	eigen_l[3][3] = -_c1_rho1_Gamma*_w + _rho1;
	eigen_l[3][4] = _c1_rho1_Gamma;

	eigen_l[4][0] = 0.5*_c1_rho1_Gamma*q2 + _w*_rho1;
	eigen_l[4][1] = -_c1_rho1_Gamma*_u;
	eigen_l[4][2] = -_c1_rho1_Gamma*_v;
	eigen_l[4][3] = -_c1_rho1_Gamma*_w - _rho1;
	eigen_l[4][4] = _c1_rho1_Gamma;

	//right eigen vectors
	eigen_r[0][0] = 0.0;
	eigen_r[0][1] = 0.0;
	eigen_r[0][2] = 1.0;
	eigen_r[0][3] = _c1_rho;
	eigen_r[0][4] = _c1_rho;
	
	eigen_r[1][0] = 0.0;
	eigen_r[1][1] = - _rho;
	eigen_r[1][2] = _u;
	eigen_r[1][3] = _c1_rho*_u;
	eigen_r[1][4] = _c1_rho*_u;

	eigen_r[2][0] = _rho;
	eigen_r[2][1] = 0.0;
	eigen_r[2][2] = _v;
	eigen_r[2][3] = _c1_rho*_v;
	eigen_r[2][4] = _c1_rho*_v;

	eigen_r[3][0] = 0.0;
	eigen_r[3][1] = 0.0;
	eigen_r[3][2] = _w;
	eigen_r[3][3] = _c1_rho*(_w + _c);
	eigen_r[3][4] = _c1_rho*(_w - _c);

	eigen_r[4][0] = _rho*_v;
	eigen_r[4][1] = -_rho*_u;
	eigen_r[4][2] = 0.5*q2;
	eigen_r[4][3] = _c1_rho*(_H + _w*_c);
	eigen_r[4][4] = _c1_rho*(_H - _w*_c);

	eigen_value[0] = fabs(_w);
	eigen_value[1] = eigen_value[0];
	eigen_value[2] = eigen_value[0];
	eigen_value[3] = fabs(_w + _c);
	eigen_value[4] = fabs(_w - _c);
	ave_rho1 = _rho1;
}
//-------------------------------------------------------------------------------------------------
//		get global maximum eigen values x direction
//-------------------------------------------------------------------------------------------------
void BaseFluid::GetMaxEigen_x()
{
	for(int n=0; n<5; n++) {
		eigen_max[n] = 0.0; //for maximum eigen values
	}

	//Eigen values
	for(int i=0; i<=l+8; i++) 
		for(int j=0; j<=d+8; j++) 
			for(int k=0; k<=width+8; k++) {
				//get the largest velocity  of the domain
				double uu = u[i][j][k];
				double uuPc = u[i][j][k] + c[i][j][k];
				double uuMc = u[i][j][k] - c[i][j][k];
				
				//local eigen values
				eigen_local[0][i][j][k] = uu;
				eigen_local[1][i][j][k] = uu;
				eigen_local[2][i][j][k] = uu;
				eigen_local[3][i][j][k] = uuPc;
				eigen_local[4][i][j][k] = uuMc;

				//global eigen values
				eigen_max[0] = AMAX1(eigen_max[0], fabs(uu));
				eigen_max[1] = eigen_max[0];
				eigen_max[2] = eigen_max[0];
				eigen_max[3] = AMAX1(eigen_max[3], fabs(uuPc));
				eigen_max[4] = AMAX1(eigen_max[4], fabs(uuMc));
			}
}
//-------------------------------------------------------------------------------------------------
//		get global maximum eigen values y direction
//-------------------------------------------------------------------------------------------------
void BaseFluid::GetMaxEigen_y()
{
	for(int n=0; n<5; n++) {
		eigen_max[n] = 0.0; //for maximum eigen values
	}
	//For local lax-friedrichs
	for(int i=0; i<=l+8; i++) 
		for(int j=0; j<=d+8; j++) 
			for(int k=0; k<=width+8; k++) {
				//get the largest velocity  of the domain
				double uu = v[i][j][k];
				double uuPc = v[i][j][k] + c[i][j][k];
				double uuMc = v[i][j][k] - c[i][j][k];
				
				//local eigen values
				eigen_local[0][i][j][k] = uu;
				eigen_local[1][i][j][k] = uu;
				eigen_local[2][i][j][k] = uu;
				eigen_local[3][i][j][k] = uuPc;
				eigen_local[4][i][j][k] = uuMc;

				//global eigen values
				eigen_max[0] = AMAX1(eigen_max[0], fabs(uu));
				eigen_max[1] = eigen_max[0];
				eigen_max[2] = eigen_max[0];
				eigen_max[3] = AMAX1(eigen_max[3], fabs(uuPc));
				eigen_max[4] = AMAX1(eigen_max[4], fabs(uuMc));
			}

}
//-------------------------------------------------------------------------------------------------
//		get global maximum eigen values z direction
//-------------------------------------------------------------------------------------------------
void BaseFluid::GetMaxEigen_z()
{
	for(int n=0; n<5; n++) {
		eigen_max[n] = 0.0; //for maximum eigen values
	}
	//For local lax-friedrichs
	for(int i=0; i<=l+8; i++) 
		for(int j=0; j<=d+8; j++) 
			for(int k=0; k<=width+8; k++) {
				//get the largest velocity  of the domain
				double uu = w[i][j][k];
				double uuPc = w[i][j][k] + c[i][j][k];
				double uuMc = w[i][j][k] - c[i][j][k];
				
				//local eigen values
				eigen_local[0][i][j][k] = uu;
				eigen_local[1][i][j][k] = uu;
				eigen_local[2][i][j][k] = uu;
				eigen_local[3][i][j][k] = uuPc;
				eigen_local[4][i][j][k] = uuMc;

				//global eigen values
				eigen_max[0] = AMAX1(eigen_max[0], fabs(uu));
				eigen_max[1] = eigen_max[0];
				eigen_max[2] = eigen_max[0];
				eigen_max[3] = AMAX1(eigen_max[3], fabs(uuPc));
				eigen_max[4] = AMAX1(eigen_max[4], fabs(uuMc));
			}

}
//-------------------------------------------------------------------------------------------------
//					Obtain state at a grid point
//-------------------------------------------------------------------------------------------------
void BaseFluid::GetStates(int i, int j, int k, double UI[Emax][Xmax][Ymax][Zmax])
{
	rho[i][j][k]	=	UI[0][i][j][k];
	double rho1		=	1.0/rho[i][j][k];
	u[i][j][k]		=	UI[1][i][j][k]*rho1;
	v[i][j][k]		=	UI[2][i][j][k]*rho1;
	w[i][j][k]		=	UI[3][i][j][k]*rho1;
	double rho_e	=	UI[4][i][j][k] - 0.5*rho[i][j][k]*(u[i][j][k]*u[i][j][k] 
						+ v[i][j][k]*v[i][j][k] + w[i][j][k]*w[i][j][k]);
	p[i][j][k]		=	_Gamma*rho_e;
	H[i][j][k]		=	(UI[4][i][j][k] + p[i][j][k])*rho1;
	c[i][j][k]		=	sqrt(Gamma*p[i][j][k]*rho1);
	e[i][j][k]		=	rho_e*rho1;
}
//-------------------------------------------------------------------------------------------------
//					Obtain conservatives at a grid point
//-------------------------------------------------------------------------------------------------
void BaseFluid::GetU(int i, int j, int k, double UI[Emax][Xmax][Ymax][Zmax])
{
	UI[0][i][j][k] = rho[i][j][k];
	UI[1][i][j][k] = rho[i][j][k]*u[i][j][k];
	UI[2][i][j][k] = rho[i][j][k]*v[i][j][k];
	UI[3][i][j][k] = rho[i][j][k]*w[i][j][k];
	UI[4][i][j][k] = _Gamma1*p[i][j][k] + 0.5*rho[i][j][k]*(u[i][j][k]*u[i][j][k] 
				   + v[i][j][k]*v[i][j][k] + w[i][j][k]*w[i][j][k]);
}
//-------------------------------------------------------------------------------------------------
//					time-step without visocus effects
//-------------------------------------------------------------------------------------------------
double BaseFluid::ViscousHeat_dt_0(double cfl_dt, double rho_min)
{
	return cfl_dt;
}
//-------------------------------------------------------------------------------------------------
//					time-step with visocus effects (not finish!)
//-------------------------------------------------------------------------------------------------
double BaseFluid::ViscousHeat_dt_1(double cfl_dt, double rho_min)
{
	double nu_alpha = AMAX1(1.0, 1.0/Pr)*mu1/rho_min;
	return AMIN1(cfl_dt, 0.5/(nu_alpha/dx/dx + nu_alpha/dy/dy + nu_alpha/dz/dz));
}
//-------------------------------------------------------------------------------------------------
//					Obtain fluxes at a grid point
//-------------------------------------------------------------------------------------------------
void BaseFluid::GetFxyz(int i, int j, int k, double UI[Emax][Xmax][Ymax][Zmax])
{
	F_x[0][i][j][k] = UI[1][i][j][k];
	F_x[1][i][j][k] = UI[1][i][j][k]*u[i][j][k] + p[i][j][k];
	F_x[2][i][j][k] = UI[1][i][j][k]*v[i][j][k];
	F_x[3][i][j][k] = UI[1][i][j][k]*w[i][j][k];
	F_x[4][i][j][k] = (UI[4][i][j][k] + p[i][j][k])*u[i][j][k];

	F_y[0][i][j][k] = UI[2][i][j][k];
	F_y[1][i][j][k] = UI[2][i][j][k]*u[i][j][k];
	F_y[2][i][j][k] = UI[2][i][j][k]*v[i][j][k] + p[i][j][k];
	F_y[3][i][j][k] = UI[2][i][j][k]*w[i][j][k];
	F_y[4][i][j][k] = (UI[4][i][j][k] + p[i][j][k])*v[i][j][k];

	F_z[0][i][j][k] = UI[3][i][j][k];
	F_z[1][i][j][k] = UI[3][i][j][k]*u[i][j][k];
	F_z[2][i][j][k] = UI[3][i][j][k]*v[i][j][k];
	F_z[3][i][j][k] = UI[3][i][j][k]*w[i][j][k] + p[i][j][k];
	F_z[4][i][j][k] = (UI[4][i][j][k] + p[i][j][k])*w[i][j][k];
}
//-------------------------------------------------------------------------------------------------
//					Initialize a empty fluid domain
//-------------------------------------------------------------------------------------------------
void SingleFluid::FluidInitializer(Initialization &Ini)
{
	char Key_word[100];
	
	//get domain dimension
	l = Ini.l; d = Ini.d; width = Ini.width; 
	lc = Ini.lc;
	dx = Ini.dx; dy = Ini.dy; dz = Ini.dz;

	//time start for computing performance
	time_start = clock();
	exec_time = time_start;

	//interpolation method
	switch(interpolation_method) {
	
		case 1: 
	   	CurrentInterpolation_P = weno5old_P;
	   	CurrentInterpolation_M = weno5old_M;
		CurrentLinearInterpolation_du = du_upwind5;
		CurrentLinearInterpolation_f2 = f2_upwind5;
		hybrid_epsilon = 1.0e2*dx*dx*dx/lc/lc/lc;
		stencil_P = 2; stencil_size = 6;
		break;

		case 2:
		CurrentInterpolation_P = weno5LP_P;
		CurrentInterpolation_M = weno5LP_M;
		CurrentLinearInterpolation_du = du_upwind5;
		CurrentLinearInterpolation_f2 = f2_upwind5;
		hybrid_epsilon = 1.0e2*dx*dx*dx / lc / lc / lc;
		stencil_P = 2; stencil_size = 6;
		break;

		case 3:
	   	CurrentInterpolation_P = wenoZ_P;
	   	CurrentInterpolation_M = wenoZ_M;
		stencil_P = 2; stencil_size = 6;
		break;

		case 4:
	   	CurrentInterpolation_P = weno_cu6_P;
	   	CurrentInterpolation_M = weno_cu6_M;
		CurrentLinearInterpolation_du = du_optimal5;
		CurrentLinearInterpolation_f2 = f2_optimal5; 
		hybrid_epsilon = 1.0e-2*dx*dx*dx / lc / lc / lc;
		stencil_P = 2; stencil_size = 6;
		break;

		case 5:
	   	CurrentInterpolation_P = weno7_P;
	   	CurrentInterpolation_M = weno7_M;
		CurrentLinearInterpolation_du = du_upwind7;
		CurrentLinearInterpolation_f2 = f2_upwind7;
		hybrid_epsilon = 1.0e-6;
		stencil_P = 3; stencil_size = 8;
		break;

		case 6:
	   	CurrentInterpolation_P = weno8_P;
	   	CurrentInterpolation_M = weno8_M;
		stencil_P = 3; stencil_size = 8;
		break;
	}

	//Lax-Friedrichs flux type
	//1: global LF, 2: local LF, 3:Roe
	switch(flux_method) {
		case 1: 
		Roe_type = 0.0; LLF_type = 0.0; GLF_type = 1.0; 
		break;

		case 2:
		Roe_type = 0.0; LLF_type = 1.0; GLF_type = 0.0; 
		break;

		case 3:
		Roe_type = 1.0; LLF_type = 0.0; GLF_type = 0.0; 
		break;
	}


	Ini.InitialCondition(this,rho, p, u, v, w);
	_Gamma = Gamma - 1.0;
	_Gamma1 = 1.0 / _Gamma;

	//input data form a restart file if needed
	if(Ini.restart_file == 1) {
		ifstream fin("singlefluid.rst");
		while(!fin.eof()) {
		
			fin>>Key_word;
			if(!strcmp(Key_word,"*conservatives")){
				
				//read state of all the grid points
				for(int i=0; i<=l+8; i++) 
					for(int j=0; j<=d+8; j++) 
						for(int k=0; k<=width+8; k++){
							for(int n=0; n<5; n++)	fin>>U[n][i][j][k];
							GetStates(i, j, k, U);
				}
				break;
			}
		}
		fin.close();
	}

	//input data form initial data for homogeneous decay turbulence
	if(Ini.restart_file == 2) {
		ifstream fin("decay_turbulence.rst");
		
		int itt = 0;
		//read state of all the grid points
		for(int n=0; n<5; n++)
			for(int k=3; k<=width+5; k++)
				for(int j=3; j<=d+5; j++)
					for(int i=3; i<=l+5; i++){
						fin>>U[n][i][j][k];
						itt ++;
		}
		fin.close();
		//boudnary condition
		Boundary(U);
		//update state
		for(int i=0; i<=l+8; i++) 
				for(int j=0; j<=d+8; j++) 
					for(int k=0; k<=width+8; k++){
						GetStates(i, j, k, U);
				}
	}

	//initialize U, States, Fluxes
	for(int i=0; i<=l+8; i++) {
		for(int j=0; j<=d+8; j++) 
			for(int k=0; k<=width+8; k++){
				GetU(i, j, k, U);
				GetStates(i, j, k, U);
				GetFxyz(i, j, k, U);
				//hybrid indicator
				hybrid_indcator[0][i][j][k] = 0.0;
				hybrid_indcator[1][i][j][k] = 0.0;
				//interval values
				for(int n=0; n<5; n++) {
					LU[n][i][j][k] = 0.0; //incremental of one time step
					U1[n][i][j][k] = U[n][i][j][k]; //intermediate conwervatives
					F_x_wall[n][i][j][k] = 0.0; //numerical flux in x direction
					F_y_wall[n][i][j][k] = 0.0; //numerical flux in y direction
					F_z_wall[n][i][j][k] = 0.0; //numerical flux in z direction
			}
		}
	}
}
//-------------------------------------------------------------------------------------------------
//					Get dt for the domain
//-------------------------------------------------------------------------------------------------
double SingleFluid::Get_dt()
{
	double uc_max=0.0, vc_max=0.0, wc_max=0.0, rho_min=1.0e13;
	for(int i=4; i<=l+4; i++) {
		for(int j=4; j<=d+4; j++) 
			for(int k=4; k<=width+4; k++){
			uc_max = AMAX1(uc_max, fabs(u[i][j][k]) + c[i][j][k]);
			vc_max = AMAX1(vc_max, fabs(v[i][j][k]) + c[i][j][k]);
			wc_max = AMAX1(wc_max, fabs(w[i][j][k]) + c[i][j][k]);
			rho_min = AMIN1(rho_min, rho[i][j][k]);
		}
	}
	double ddt = 1.0/(uc_max/dx + vc_max/dy + wc_max/dz); 
	dt = (this->*CurrentViscousHeat_dt)(cfl*ddt, rho_min);
	lambda_x = cfl/uc_max; lambda_y = cfl/vc_max; lambda_z = cfl/wc_max; 
	lambda_x0 = uc_max; lambda_y0 = vc_max; lambda_z0 = wc_max; 
	return dt;
}
//-------------------------------------------------------------------------------------------------
//	Update U form pressure, density and velocity
//-------------------------------------------------------------------------------------------------
void SingleFluid::Renew_U(double UI[Emax][Xmax][Ymax][Zmax])
{
	for(int i=0; i<=l+8; i++) {
		for(int j=0; j<=d+8; j++) 
			for(int k=0; k<=width+8; k++){
					GetU(i, j, k, UI);
		}
	}
}
//-------------------------------------------------------------------------------------------------
//		Renew states from conservative variables
//-------------------------------------------------------------------------------------------------
void  SingleFluid::Renew_state(double UI[Emax][Xmax][Ymax][Zmax])
{
	for(int i=0; i<=l+8; i++) {
		for(int j=0; j<=d+8; j++) 
			for(int k=0; k<=width+8; k++){
				GetStates(i, j, k, UI);
				GetFxyz(i, j, k, UI);
		}
	}
}
//-------------------------------------------------------------------------------------------------
//		Re_new the fluxes
//-------------------------------------------------------------------------------------------------
void SingleFluid::Renew_Fluxes(double UI[Emax][Xmax][Ymax][Zmax])
{
	double epsilon = 1.0e-13;
	//proceed at x directiom and get F-flux terms at node wall
	GetMaxEigen_x();
	time_now = clock();
	for (int j = 4; j <= d + 4; j++)
	for (int k = 4; k <= width + 4; k++)
		//get the F-flux terms at cell wall
	for (int i = 3; i <= l + 4; i++){
		CharacteristicF_x(i, j, k, UI);
	}
	exec_time_interval = clock() - time_now;
//	cout<<exec_time_interval<<"\n";
	exec_time += exec_time_interval;
		
	if(d != 0) {
		//proceed at y directiom and get G-flux terms at node wall
		GetMaxEigen_y();
		time_now = clock();
		for (int i = 4; i <= l + 4; i++)
			for(int k=4; k<=width+4; k++)
				//get the G-flux terms at cell wall
				for(int j=3; j<=d+4; j++){
					CharacteristicF_y(i, j, k, UI);
		}
		exec_time_interval = clock() - time_now;
		exec_time += exec_time_interval;
	}

	if(width != 0) {
		//proceed at z directiom and get G-flux terms at node wall
		GetMaxEigen_z();
		for(int i=4; i<=l+4; i++)
			for(int j=4; j<=d+4; j++)
				//get the G-flux terms at cell wall
				for(int k=3; k<=width+4; k++){
					CharacteristicF_z(i, j, k, UI);
		}	
	}	
	//viscous effect
	(this->*CurrentPhysicalDissipation)();
	//positivity preserving
	(this->*CurrentPositivityPreserving)(UI);

	//get LU
	for(int i=4; i<=l+4; i++)
		for(int j=4; j<=d+4; j++) 
			for(int k=4; k<=width+4; k++) 
				for(int n=0; n<5; n++) {
					LU[n][i][j][k] = (F_x_wall[n][i-1][j][k] - F_x_wall[n][i][j][k])/dx
								   + (F_y_wall[n][i][j-1][k] - F_y_wall[n][i][j][k])/dy
								   + (F_z_wall[n][i][j][k-1] - F_z_wall[n][i][j][k])/dz
								   ;
				}
}
//-------------------------------------------------------------------------------------------------
//		March one sub-time step with 3rd Runge_kutta
//		flag; 1 first, 2 second and 3 third sub-time step
//-------------------------------------------------------------------------------------------------
void SingleFluid::Step_rk3(double Time, double dt, int flag)
{
	switch(flag) {
	
	case 1:
		//the fisrt step
		run_time = Time;
		Boundary(U);
		Renew_state(U);
		Renew_Fluxes(U);
		for(int i=4; i<=l+4; i++) 
			for(int j=4; j<=d+4; j++)
				for(int k=4; k<=width+4; k++)
					for(int n=0; n<5; n++) {
						U1[n][i][j][k] = U[n][i][j][k] + dt*LU[n][i][j][k];
					}
	break;
	
	
	case 2:
		//the second step
		Boundary(U1);
		Renew_state(U1);
		Renew_Fluxes(U1);
		for(int i=4; i<=l+4; i++) 
			for(int j=4; j<=d+4; j++) 
				for(int k=4; k<=width+4; k++)
					for(int n=0; n<5; n++) {
						U1[n][i][j][k] = 0.75*U[n][i][j][k] + 0.25*U1[n][i][j][k] 
										 + 0.25*dt*LU[n][i][j][k];
				}
	break;


	case 3:
		//the third step
		Boundary(U1);
		Renew_state(U1);
		Renew_Fluxes(U1);
		for(int i=4; i<=l+4; i++) 
			for(int j=4; j<=d+4; j++) 
				for(int k=4; k<=width+4; k++)
					for(int n=0; n<5; n++) {
						U[n][i][j][k] = (U[n][i][j][k] + 2.0*U1[n][i][j][k])/3.0
										+ 2.0*dt*LU[n][i][j][k]/3.0;
				}
		break;
	}
}
