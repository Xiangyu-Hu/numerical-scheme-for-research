const double pi=3.1415926535897932384626434; //Pi
//	a. Get the maximum
double AMAX1(double a, double b);
//	b. Get the minimum
double AMIN1(double a, double b);
//	c. Get the sign of the a double
double sgn(double a);
//	d. Sign of the first value is determined by the secomax_species_number value's sign
double SIGN(double a, double b);
//	e. Get the maximum from four
double AMAX4(double a, double b, double c, double d);
//	f. Get the minimum from four
double AMIN4(double a, double b, double c, double d);
//	g. minmod function
double minmod(double a, double b);
//	h. minmod function
double median(double a, double b, double c);
//  i. Heaviside function
double Heaviside(double phi, double epsilon);
//  j. Wenland function
double Wenland(double phi, double epsilon);
//----------------------------------------------------------------------------------------
//			get the random number with uniform distribution in [0, 1]
//----------------------------------------------------------------------------------------
double Ranuls();
//----------------------------------------------------------------------------------------
//		get three random numbers y1, y2, y3 with guassian distribution 
//		with zero mean and variance one
//		from random numbers uniform distributed in [0, 1]
//----------------------------------------------------------------------------------------
void  Gaussian(double &y1, double &y2, double &y3);
//---------------------------------------------------------------------------------
//	function for calculate the nth order divided difference
//	referrence 水鸿寿：一维流体力学差分方法，国防工业出版社，第480页
//	f: the adress of value on the node j
//	n: the order of divided difference
//	delta: space steps such as dx or dy
//	return: the nth order divided difference at j(for even orders), or j+1/2(for odd oders)
//---------------------------------------------------------------------------------
double div_diff(double *f, int n, double delta);
//---------------------------------------------------------------------------------
//	curvature measurement
//---------------------------------------------------------------------------------
double curv(double u_p, double u, double u_m);
double M4X(double a, double b, double c, double d, double e, double f);
//---------------------------------------------------------------------------------
//         ENO scheme
//---------------------------------------------------------------------------------
//	function for calculating the numerical flux of scalar equations with 3th eno
//	referrence: Ronald P. Fedkiw et al., 
//	and Chi-Wang Shu et al.,  Efficient Implementation of Essential Non-osillatory Shock-
//	Capturing Schemes, II, JCP, vol 83, 32-78, 1989.
//	lambda,: the maximum wave speeds 
//	delta: space steps such as dx or dy
//	return: the numerical flux at j+1/2
//---------------------------------------------------------------------------------
//this is ENO-LLF which is used for rarefication wave
double eno_P(double *f, double delta);
double eno_M(double *f, double delta);
// the end of scalar ENO scheme
//-----------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//         WENO scheme
//---------------------------------------------------------------------------------
//	function for calculating the numerical flux of scalar equations with 3th eno
//	referrence: Ronald P. Fedkiw et al., 
//	and Chi-Wang Shu et al.,  Efficient Implementation of Essential Non-osillatory Shock-
//	Capturing Schemes, II, JCP, vol 83, 32-78, 1989.
//	lambda,: the maximum wave speeds 
//	delta: space steps such as dx or dy
//	return: the numerical flux at j+1/2
//---------------------------------------------------------------------------------
//this is WENO-LLF which is used for rarefication wave
double weno5old_P(double *f, double delta);
double weno5old_M(double *f, double delta);
double du_upwind5(double *f, double delta);
double f2_upwind5(double *f, double delta);
double du_upwind7(double *f, double delta);
double f2_upwind7(double *f, double delta);
double weno_cu6_P(double *f, double delta);
double weno_cu6_M(double *f, double delta);
double weno_cu6_M1_P(double *f, double delta);
double weno_cu6_M1_M(double *f, double delta);
double du_optimal5(double *f, double delta);
double f2_optimal5(double *f, double delta);
double wenoZ_P(double *f, double delta);
double wenoZ_M(double *f, double delta);
double weno5LP_P(double *f, double delta);
double weno5LP_M(double *f, double delta);

// the end of scalar ENO scheme
//-----------------------------------------------------------------------------------------
//Lower 3rd order WENO 
//-----------------------------------------------------------------------------------------
double weno3_P(double *f, double delta);
double weno3_M(double *f, double delta);
double weno3old_P(double *f, double delta);
double weno3old_M(double *f, double delta);
double weno2_P(double *f, double delta);
double weno2_M(double *f, double delta);
double weno4_P(double *f, double delta);
double weno4_M(double *f, double delta);
double weno6_P(double *f, double delta);
double weno6_M(double *f, double delta);
double weno7_P(double *f, double delta);
double weno7_M(double *f, double delta);
double weno8_P(double *f, double delta);
double weno8_M(double *f, double delta);
double weno3new_P(double *f, double delta);
double weno3new_M(double *f, double delta);
//-----------------------------------------------------------------------------------------
//		fifth order upwind central scheme 
//-----------------------------------------------------------------------------------------
double houc5_P(double *f, double delta);
double houc5_M(double *f, double delta);
//-----------------------------------------------------------------------------------------
//		fifth order WENO with uniform accuracy in smooth regions
//-----------------------------------------------------------------------------------------
double weno5u_P(double *f, double delta);
double weno5u_M(double *f, double delta);
//-----------------------------------------------------------------------------------------
//		sixth order WENO with uniform accuracy in smooth regions
//-----------------------------------------------------------------------------------------
double weno4u_P(double *f, double delta);
double weno4u_M(double *f, double delta);
//-----------------------------------------------------------------------------------------
//		sixth order center scheme
//-----------------------------------------------------------------------------------------
double central_6(double *f, double delta);
//-----------------------------------------------------------------------------------------
//		sixth order center scheme for derivative
//-----------------------------------------------------------------------------------------
double central_d_6(double *f, double delta);
//-----------------------------------------------------------------------------------------
//		resolvable index
//-----------------------------------------------------------------------------------------
double index_r(double *f, double delta);
//-----------------------------------------------------------------------------------------
//		fifth order WENO with uniform accuracy in smooth regions
//-----------------------------------------------------------------------------------------
double weno5u_P(double *f, double delta, double r_index);
double weno5u_M(double *f, double delta, double r_index);
