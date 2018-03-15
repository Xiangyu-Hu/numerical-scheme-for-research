//-------------------------------------------------------------------------------------------------
//							Basic global values
//-------------------------------------------------------------------------------------------------
const int Xmax = 810; // maximum number of nodes in x direction
const int Ymax = 10;  // maximum number of nodes in y direction
const int Zmax = 10;  // maximum number of nodes in Z direction
const int Emax = 5; // maximum number of equations, 0 density, 1 x_momentum, 2 y_momentum, 3 z_momentum, 4 energy
//-------------------------------------------------------------------------------------------------
//								Pre-claimer
//-------------------------------------------------------------------------------------------------
class BaseSolver; class Output;
class BaseMaterial; class BaseFluid; 
class SingleFluid; class BaseSolver;
//--------------------------------------------------------------------------------------------------
//					Member date and memebr function poiters
//--------------------------------------------------------------------------------------------------
typedef double (*Interpolation_FP)(double *f, double delta);
typedef void (BaseFluid::*BaseFluidFluxMethod_FP)(int i, int j, int k, double UI[Emax][Xmax][Ymax][Zmax]);
typedef void (BaseFluid::*BaseFluidPhysicalDissipation_FP)();
typedef void (BaseFluid::*BaseFluidPositivityPreserving_FP)(double UI[Emax][Xmax][Ymax][Zmax]);
typedef double (BaseFluid::*BaseFluidViscousHeat_dt_FP)(double cfl_dt, double rho_min);
//-------------------------------------------------------------------------------------------------
//								Initialization
//-------------------------------------------------------------------------------------------------
class Initialization {

	//for random number generator
	//variables
	long int idum, idum2, iy, ntab, iseed;
	long int *iv;
	//set the random seed
	void Ranils();
	//get the random number uniform distributed in [0, 1]
	double Ranuls();
	//get two random numbers y1, y2 with guassian distribution with zero mean and variance one
	//from random numbers uniform distributed in [0, 1]
	void Gaussian(double &y1, double &y2, double &y3);

public:
	int l, d, width; //step numbers on x, y and z direction
	int ld, lw; //step number ration in y and z direction: should be dividable for l
	double length, dx, dy, dydx, dz, dzdx, lc; //space steps on x and y direction
	double start_time, end_time, output_period;
	int restart_file; //1: input from data file, 0: not
	int case_number, problem_number;

	Initialization(); //constructor, output the universal parameters to screen too
	void InitialCondition(BaseFluid *ActiveFluid, double rho[Xmax][Ymax][Zmax], double p[Xmax][Ymax][Zmax],
	double u[Xmax][Ymax][Zmax], double v[Xmax][Ymax][Zmax], double w[Xmax][Ymax][Zmax]);
};
//-------------------------------------------------------------------------------------------------
//							Fluid Compressible
//-------------------------------------------------------------------------------------------------
class BaseFluid{
	friend class Initialization;

protected:
	int    l, d, width; //size of the computational domain
	double dx, dy, dz, lc; //spatial step
	double cfl, dt, lambda_x, lambda_y, lambda_z, lambda_x0, lambda_y0, lambda_z0; //dt/dx, dt/dy and dt/dy;
	int xBl, xBr, yBu, yBd, zBn, zBf; //boundary condition indicator
	int interpolation_method, flux_method;//flux method: 1 global Lax-Friedrischs, 2 local Lax-Friedrischs, 3 Roe
	double g_x, g_y, g_z; //gravity
	double mu1, mu2, Pr, rho_alpha; //first and second viscosity and Prandtl number
	double run_time; //runtime

	//constant for ideal gas EOS
	double Gamma; //specific heat ratio 
	double _Gamma; //Gamma -1
	double _Gamma1; //1 / (Gamma -1)

	int stencil_P, stencil_size; //stencil span used for weno interpolation
	double Roe_type, LLF_type, GLF_type;//artificial_viscosity types
	//artificial_viscosity for global and local Lax Friedrischs types
	double eigen_max[Emax];
	double eigen_local[Emax][Xmax][Ymax][Zmax]; 
	//sound speed, pressure derivatives on conservative variables
	double c[Xmax][Ymax][Zmax], H[Xmax][Ymax][Zmax], e[Xmax][Ymax][Zmax];
	//Fluxes calculated on cell averages		
	double F_x[Emax][Xmax][Ymax][Zmax], F_y[Emax][Xmax][Ymax][Zmax], F_z[Emax][Xmax][Ymax][Zmax]; 
	double F_x_wall[Emax][Xmax][Ymax][Zmax], F_y_wall[Emax][Xmax][Ymax][Zmax], F_z_wall[Emax][Xmax][Ymax][Zmax]; 
	double hybrid_epsilon, hybrid_indcator[2][Xmax][Ymax][Zmax];

	void Boundary(double UI[Emax][Xmax][Ymax][Zmax]); //boundary condition
	void GetMaxEigen(); //get local and global maximum eigen values
	void GetMaxEigen_x(); //get local and global maximum eigen values
	void GetMaxEigen_y(); //get local and global maximum eigen values
	void GetMaxEigen_z(); //get local and global maximum eigen values
	//get Roe average and eigenvectors
	void RoeAverage_x(int i, int j, int k, double eigen_l[Emax][Emax], 
		double eigen_r[Emax][Emax], double eigen_value[Emax], double &ave_rho1); 
	void RoeAverage_y(int i, int j, int k, double eigen_l[Emax][Emax], 
		double eigen_r[Emax][Emax], double eigen_value[Emax], double &ave_rho1); 
	void RoeAverage_z(int i, int j, int k, double eigen_l[Emax][Emax], 
		double eigen_r[Emax][Emax], double eigen_value[Emax], double &ave_rho1); 
	void GetStates(int i, int j, int k, double UI[Emax][Xmax][Ymax][Zmax]); //get states from conservative variables
	void GetU(int i, int j, int k, double UI[Emax][Xmax][Ymax][Zmax]); //get conservative variables from states
	void GetFxyz(int i, int j, int k, double UI[Emax][Xmax][Ymax][Zmax]); //get cell center fluxes
	void CharacteristicF_x(int i, int j, int k, double UI[Emax][Xmax][Ymax][Zmax]); //get cell face fluxes in x direction
	void CharacteristicF_y(int i, int j, int k, double UI[Emax][Xmax][Ymax][Zmax]);//get cell face fluxes in y direction
	void CharacteristicF_z(int i, int j, int k, double UI[Emax][Xmax][Ymax][Zmax]);//get cell face fluxes in y direction
	
	double ViscousHeat_dt_0(double cfl_dt, double rho_min); //time-step without visocus effects
	double ViscousHeat_dt_1(double cfl_dt, double rho_min); //time-step for viscous flow
	void PhysicalDissipation_0(); //no visocus effects
	void PhysicalDissipation_1(); //with visocus effects
	void PositivityPreserving_0(double UI[Emax][Xmax][Ymax][Zmax]);
	void PositivityPreserving_1(double UI[Emax][Xmax][Ymax][Zmax]);

	//function pointers
	Interpolation_FP CurrentInterpolation_P, CurrentInterpolation_M;
	Interpolation_FP CurrentLinearInterpolation_du, CurrentLinearInterpolation_f2;
	BaseFluidViscousHeat_dt_FP CurrentViscousHeat_dt;
	BaseFluidPhysicalDissipation_FP CurrentPhysicalDissipation;
	BaseFluidPositivityPreserving_FP CurrentPositivityPreserving;

	virtual void Renew_state(double UI[Emax][Xmax][Ymax][Zmax]) = 0;	//renew states from conservative variables		
	virtual void Renew_U(double UI[Emax][Xmax][Ymax][Zmax]) = 0; //update the U from states
	
public:
	
	char Fluid_name[100],  mtl_name[100]; //names of the fluid and its material
	//Consevative variable, states for public usage
	double U[Emax][Xmax][Ymax][Zmax], rho[Xmax][Ymax][Zmax], p[Xmax][Ymax][Zmax];
	double u[Xmax][Ymax][Zmax], v[Xmax][Ymax][Zmax], w[Xmax][Ymax][Zmax];
	double U1[Emax][Xmax][Ymax][Zmax], LU[Emax][Xmax][Ymax][Zmax];		
	//measure execute time
	unsigned int exec_time, exec_time_interval, time_now, time_start;

	BaseFluid(); //constructor
	//for output 
	void OutputStates(int i, int j, int k, double &denisty, double &pressure, 
							 double &uu, double &vv,  double &ww, double &Hybrid_factor); 

	virtual double Get_dt() = 0; //get dt for fluid integeral
	virtual void Step_rk3(double Time, double dt, int flag) = 0; //time integeral with 3rd Runge-Kutta	
};
//----------------------------------------------------------------------------------------------
//		For single fluid applications
//----------------------------------------------------------------------------------------------
class SingleFluid : public BaseFluid
{
	virtual void Renew_U(double UI[Emax][Xmax][Ymax][Zmax]); //update the U from states
	virtual void Renew_state(double UI[Emax][Xmax][Ymax][Zmax]);	//renew states		
	virtual void Renew_Fluxes(double UI[Emax][Xmax][Ymax][Zmax]); //update the fluxes

public:
	SingleFluid(){}; //constructor
	//initialize for a empty fluid
	void FluidInitializer(Initialization &Ini); 

	virtual double Get_dt(); //get time step for fluid integeral
	virtual void Step_rk3(double Time, double dt, int flag); //with separated sub-time Runge-Kutta steps
};
//-------------------------------------------------------------------------------------------------
//										Out_put
//-------------------------------------------------------------------------------------------------
class Output {

	double dx, dy, dz; //spatial step
	int l, d, width; //size of the computational domain
	double Tmass10, Tmass20, Tmass0, Tenergy10, Tenergy20,Tenergy0;

public:

	Output(double Time, Initialization &Ini); //constructor
	void OutStates(BaseSolver *CurrentSolver, double Time, Initialization &Ini);
	//error for case 21 problem 8
	void OutError218(BaseSolver *CurrentSolver, double Time, Initialization &Ini); 
};
//-------------------------------------------------------------------------------------------------
//							Base solver class
//-------------------------------------------------------------------------------------------------
class BaseSolver{

protected:

	int Ite; //number of iteration
	double Integeral_time; //current integeral time
	double dx, dy, dz, dt0; //spatial and time step, 
	int l, d, width; //size of the computational domain

	int CheckPoints[2]; //number of time steps for stepwise output

	//stepwise outputs
	virtual void Out_reinitiation(double Time, Initialization &Ini) = 0;

public:
	BaseSolver(); //constructor

	//type of the solver
	virtual void Time_integeral(Initialization &Ini, double &Time, double output_period) = 0;

	//get states
	virtual void OutputStates(int i, int j, int k, double &rho, double &p, 
					double &u, double &v, double &w, double &Hybrid_factor) = 0;
};
//-------------------------------------------------------------------------------------------------
//							Single-fluid solver
//-------------------------------------------------------------------------------------------------
class SingleFluidSolver : public BaseSolver
{
	
	double TotalKE0, THybrid_factor0, Max_U0, TotalVolume;
	SingleFluid fluid;

protected:
	SingleFluid *CurrentFluid;
	//stepwise outputs
	virtual void Out_GlobalValues(double Time, Initialization &Ini);

public:

	SingleFluidSolver(){}; //empty constructor
	SingleFluidSolver(Initialization &Ini); //constructor
	virtual void Out_reinitiation(double Time, Initialization &Ini);
	virtual void Time_integeral(Initialization &Ini, double &Time, double output_period); //for fluid-fluid interaction

	//get states
	virtual void OutputStates(int i, int j, int k, double &rho, double &p, 
					double &u, double &v, double &w, double &Hybrid_factor);

	//stepwise outputs
	virtual void InitializeGlobalValues(double Time, Initialization &Ini);

};
