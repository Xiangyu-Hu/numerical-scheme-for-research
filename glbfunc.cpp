#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string.h>
#include <cmath>
#include "glbfunc.h"

using namespace std;

//	a. Get the maximum
double AMAX1(double a, double b){
	if(a > b) return a;
	return b;
}
//	b. Get the minimum
double AMIN1(double a, double b){
	if(a < b) return a;
	return b;
}
//	c. Get the sign of the a double
double sgn(double a)
{
	if(a > 0.0) return 1.0;
	return -1.0;
}
//	d. Sign of the first value is determined by the secomax_species_number value's sign
double SIGN(double a, double b){
	if(b > 0.0) return fabs(a);
	return -fabs(a);
}
// e. Get the maximum from four
double AMAX4(double a, double b, double c, double d){
	return AMAX1(AMAX1(a, b), AMAX1(c, d));
}
//	f. Get the minimum from four
double AMIN4(double a, double b, double c, double d){
	return AMIN1(AMIN1(a, b), AMIN1(c, d));
}
//	g. minmod function
double minmod(double a, double b){
	return 0.5*(sgn(a) + sgn(b))*AMIN1(fabs(a), fabs(b));
}
//	h. minmod function
double median(double a, double b, double c){
	return a + minmod(b - a, c - a);
}
//  i. Heaviside function
double Heaviside(double phi, double epsilon){
	if(phi >= epsilon) return 1.0;
	if(phi <= -epsilon) return 0.0;
	return 0.5*(1.0 + phi/epsilon + sin(3.14159265*phi/epsilon)/3.14159265);
}
//  j. Wenland function
double Wenland(double phi, double epsilon) {
	double t = phi / epsilon;
	if ( t > 1.0) return 0.0;
	if (t < 0.0) return 1.0;
	double t3 = (1.0 - t)*(1.0 - t)*(1.0 - t);
	return t3*t3*(1.0 + 6.0*t + 35.0*t*t / 3.0);
}

//-------------------------------------------------------------------------------------------------
//	function for calculate the nth order divided difference
//	referrence Ë®ºèÊÙ£ºÒ»Î¬Á÷ÌåÁ¦Ñ§²î·Ö·½·¨£¬¹ú·À¹¤Òµ³ö°æÉç£¬µÚ480Ò³
//	f: the adress of value on the node j
//	n: the order of divided difference
//	delta: space steps such as dx or dy
//	return: the nth order divided difference at j(for even orders), or j+1/2(for odd oders)
//-------------------------------------------------------------------------------------------------
double div_diff(double *f, int n, double delta)
{
	int r;
	if(n >0){
		r = n % 2; //	get the remainder
		return (div_diff(f+r, n-1, delta) - div_diff(f+r-1, n-1, delta))/(n*delta);
	}
	return *f;
}
//---------------------------------------------------------------------------------
//	curvature measurement
//---------------------------------------------------------------------------------
double curv(double u_p, double u, double u_m) {

	return u_p - 2.0*u + u_m;
}
double M4X(double a, double b, double c, double d, double e, double f){
	return 0.5*(sgn(a) + sgn(b))
			*fabs(0.5*(sgn(a) + sgn(c))*0.5*(sgn(a) + sgn(d))*0.5*(sgn(a) + sgn(e))*0.5*(sgn(a) + sgn(f)))
			*AMIN4(AMIN1(fabs(a), fabs(b)), AMIN1(fabs(c), fabs(d)), fabs(e), fabs(f))/6.0;
}
//-------------------------------------------------------------------------------------------------
//	function for calculating the numerical flux of scalar equations with 3th eno
//	referrence: Ronald P. Fedkiw et al., The Penultimate Scheme for System of Conservation
//	laws: Finite Difference ENO with Marquina's Flux Splitting, http://www.math.ucla.edu/
//	applied/cam/index.html 
//	and Chi-Wang Shu et al.,  Efficient Implementation of Essential Non-osillatory Shock-
//	Capturing Schemes, II, JCP, vol 83, 32-78, 1989.
//	uR,uL: the adress of the  scalar conservative terms value at j+1/2 evaluated from the right and the left
//	f: the adress of value on the node j 
//	lambdaR, lambdaL: the sound speeds which evaluated from the right and the left
//	delta: space steps such as dx or dy
//	return: the numerical flux at j+1/2
//-------------------------------------------------------------------------------------------------
//this is ENO-LLF
double eno_P(double *f, double delta)
{
	int k;
	double dQ1, dQ2, dQ2L, dQ3, dQ3L, ex;
	double AdQ2, AdQ2L, AdQ3, AdQ3L;

	k = 0;
	dQ1 = div_diff(f+k, 0, delta);
	
	ex = 0.5*(-2.0*k + 1.0)*delta;
	dQ2 = div_diff(f+k, 1, delta); 
	dQ2L = div_diff(f+k-1, 1, delta);
	AdQ2 = dQ2*dQ2; AdQ2L = dQ2L*dQ2L;
	AdQ2 = AdQ2*AdQ2; AdQ2L = AdQ2L*AdQ2L; 
	if(AdQ2 > AdQ2L) k = k - 1; 
	dQ2 = (dQ2*AdQ2L + dQ2L*AdQ2)/(AdQ2 + AdQ2L + 1.0e-30);
	dQ2 = dQ2*ex;

	ex = (3.0*k*k - 1.0)*delta*delta/3.0;
	dQ3 = div_diff(f+k+1, 2, delta); 
	dQ3L =div_diff(f+k, 2, delta);
	AdQ3 = dQ3*dQ3*dQ3*dQ3; AdQ3L = dQ3L*dQ3L*dQ3L*dQ3L; 
	AdQ3 = AdQ3*AdQ3; AdQ3L = AdQ3L*AdQ3L; 
	dQ3 = (dQ3*AdQ3L + dQ3L*AdQ3)/(AdQ3 + AdQ3L + 1.0e-30);
	dQ3 = dQ3*ex;
	
	return dQ1 + dQ2 + dQ3;	
}
//
double eno_M(double *f, double delta)
{
	int k;
	double dQ1, dQ2, dQ2L, dQ3, dQ3L, ex;
	double AdQ2, AdQ2L, AdQ3, AdQ3L;

	k = 1;
	dQ1 = div_diff(f+k, 0, delta);
	
	ex = 0.5*(-2.0*k + 1.0)*delta;
	dQ2 = div_diff(f+k, 1, delta); 
	dQ2L = div_diff(f+k-1, 1, delta);
	AdQ2 = dQ2*dQ2; AdQ2L = dQ2L*dQ2L;
	AdQ2 = AdQ2*AdQ2; AdQ2L = AdQ2L*AdQ2L; 
	if(AdQ2 > AdQ2L) k = k - 1; 
	dQ2 = (dQ2*AdQ2L + dQ2L*AdQ2)/(AdQ2 + AdQ2L + 1.0e-30);
	dQ2 = dQ2*ex;

	ex = (3.0*k*k - 1.0)*delta*delta/3.0;
	dQ3 = div_diff(f+k+1, 2, delta); 
	dQ3L =div_diff(f+k, 2, delta);
	AdQ3 = dQ3*dQ3*dQ3*dQ3; AdQ3L = dQ3L*dQ3L*dQ3L*dQ3L; 
	AdQ3 = AdQ3*AdQ3; AdQ3L = AdQ3L*AdQ3L; 
	dQ3 = (dQ3*AdQ3L + dQ3L*AdQ3)/(AdQ3 + AdQ3L + 1.0e-30);
	dQ3 = dQ3*ex;
	
	return dQ1 + dQ2 + dQ3;	
}
// the end of scalar ENO scheme
//-------------------------------------------------------------------------------------------------
//       the 5th WENO Scheme
//  reference Appendix of paper: 
//A non-osillatory Eulerian Approch to interfaces in Multimaterial flows
//(Ghost Fluid Method) JCP vol 152, p457-492
//-------------------------------------------------------------------------------------------------
double weno5old_P(double *f, double delta)
{

	int k;
	double v1, v2, v3, v4, v5;
	double a1, a2, a3, w1, w2, w3;

	//assign value to v1, v2,...
	k = 0;
	v1 = *(f + k - 2);
	v2 = *(f + k - 1);
	v3 = *(f + k);
	v4 = *(f + k + 1); 
	v5 = *(f + k + 2);

	//smoothness indicator
	double s1 = 13.0*(v1 - 2.0*v2 + v3)*(v1 - 2.0*v2 + v3) 
	   + 3.0*(v1 - 4.0*v2 + 3.0*v3)*(v1 - 4.0*v2 + 3.0*v3);
	double s2 = 13.0*(v2 - 2.0*v3 + v4)*(v2 - 2.0*v3 + v4) 
	   + 3.0*(v2 - v4)*(v2 - v4);
	double s3 = 13.0*(v3 - 2.0*v4 + v5)*(v3 - 2.0*v4 + v5) 
	   + 3.0*(3.0*v3 - 4.0*v4 + v5)*(3.0*v3 - 4.0*v4 + v5);

	//weights
	a1 = 0.1/((1.0e-6 + s1)*(1.0e-6 + s1));
	a2 = 0.6/((1.0e-6 + s2)*(1.0e-6 + s2));
	a3 = 0.3/((1.0e-6 + s3)*(1.0e-6 + s3));
	double tw1 = 1.0/(a1 + a2 +a3); 
	w1 = a1*tw1;
	w2 = a2*tw1;
	w3 = a3*tw1;

	//return weighted average
	return  w1*(2.0*v1 - 7.0*v2 + 11.0*v3)/6.0
		  + w2*(-v2 + 5.0*v3 + 2.0*v4)/6.0
		  + w3*(2.0*v3 + 5.0*v4 - v5)/6.0;

}
double weno5old_M(double *f, double delta)
{

	int k;
	double v1, v2, v3, v4, v5;
	double a1, a2, a3, w1, w2, w3;

	//assign value to v1, v2,...
	k = 1;
	v1 = *(f + k + 2);
	v2 = *(f + k + 1);
	v3 = *(f + k);
	v4 = *(f + k - 1); 
	v5 = *(f + k - 2);

	//smoothness indicator
	double s1 = 13.0*(v1 - 2.0*v2 + v3)*(v1 - 2.0*v2 + v3) 
	   + 3.0*(v1 - 4.0*v2 + 3.0*v3)*(v1 - 4.0*v2 + 3.0*v3);
	double s2 = 13.0*(v2 - 2.0*v3 + v4)*(v2 - 2.0*v3 + v4) 
	   + 3.0*(v2 - v4)*(v2 - v4);
	double s3 = 13.0*(v3 - 2.0*v4 + v5)*(v3 - 2.0*v4 + v5) 
	   + 3.0*(3.0*v3 - 4.0*v4 + v5)*(3.0*v3 - 4.0*v4 + v5);

	//weights
	a1 = 0.1/((1.0e-6 + s1)*(1.0e-6 + s1));
	a2 = 0.6/((1.0e-6 + s2)*(1.0e-6 + s2));
	a3 = 0.3/((1.0e-6 + s3)*(1.0e-6 + s3));
	double tw1 = 1.0/(a1 + a2 +a3); 
	w1 = a1*tw1;
	w2 = a2*tw1;
	w3 = a3*tw1;

	//return weighted average
	return  w1*(2.0*v1 - 7.0*v2 + 11.0*v3)/6.0
		  + w2*(-v2 + 5.0*v3 + 2.0*v4)/6.0
		  + w3*(2.0*v3 + 5.0*v4 - v5)/6.0;
}
//-------------------------------------------------------------------------------------------------
//		linear scheme for hybrid method
//-------------------------------------------------------------------------------------------------
double du_upwind5(double *f, double delta)
{
	double v1 = *(f - 2);
	double v2 = *(f - 1);
	double v3 = *f;
	double v4 = *(f + 1); 
	double v5 = *(f + 2);
	double v6 = *(f + 3);
	return (v1 - 5.0*v2 + 10.0*v3 - 10.0*v4 + 5.0*v5 - v6)/60.0;
}
double f2_upwind5(double *f, double delta)
{
	double v1 = *(f - 2);
	double v2 = *(f - 1);
	double v3 = *f;
	double v4 = *(f + 1); 
	double v5 = *(f + 2);
	double v6 = *(f + 3);
	return (v1 - 8.0*v2 + 37.0*v3 + 37.0*v4 - 8.0*v5 + v6)/60.0;
}//-------------------------------------------------------------------------------------------------
//  the 6th order WENO Scheme
//  reference Paper: Hu, Wang and Adams,
//A adaptive central-upwind 6th order WENO scheme.
//JCP vol 229 (2010) 8952–8965
//-------------------------------------------------------------------------------------------------
double weno_cu6_P(double *f, double delta)
{
	//assign value to v1, v2,...
	int k = 0;
	double v1 = *(f + k - 2);
	double v2 = *(f + k - 1);
	double v3 = *(f + k);
	double v4 = *(f + k + 1); 
	double v5 = *(f + k + 2);
	double v6 = *(f + k + 3);

	//smoothness indicator
	double epsilon = 1.e-20;
	double s11 = v1 - 2.0*v2 + v3;
	double s12 = v1 - 4.0*v2 + 3.0*v3;
	double s1 = 13.0*s11*s11 + 3.0*s12*s12;
	double s21 = v2 - 2.0*v3 + v4;
	double s22 = v2 - v4;
	double s2 = 13.0*s21*s21 + 3.0*s22*s22;
	double s31 = v3 - 2.0*v4 + v5;
	double s32 = 3.0*v3 - 4.0*v4 + v5;
	double s3 = 13.0*s31*s31 + 3.0*s32*s32;
	double tau61 = (259.0*v6 - 1895.0*v5 + 6670.0*v4 - 2590.0*v3 - 2785.0*v2 + 341.0*v1) / 5760.0;
	double tau62 = -(v5 - 12.0*v4 + 22.0*v3 - 12.0*v2 + v1) / 16.0;
	double tau63 = -(7.0*v6 - 47.0*v5 + 94.0*v4 - 70.0*v3 + 11.0*v2 + 5.0*v1) / 144.0;
	double tau64 = (v5 - 4.0*v4 + 6.0*v3 - 4.0*v2 + v1) / 24.0;
	double tau65 = -(-v6 + 5.0*v5 - 10.0*v4 + 10.0*v3 - 5.0*v2 + v1) / 120.0;
	double a1a1 = 1.0, a2a2 = 13.0 / 3.0, a1a3 = 0.5, a3a3 = 3129.0 / 80.0, a2a4 = 21.0 / 5.0;
	double a1a5 = 1.0 / 8.0, a4a4 = 87617.0 / 140.0, a3a5 = 14127.0 / 224.0, a5a5 = 252337135.0 / 16128.0;
	double s6 = (tau61*tau61*a1a1 + tau62*tau62*a2a2 + tau61*tau63*a1a3 + tau63*tau63*a3a3 + tau62*tau64*a2a4
		+ tau61*tau65*a1a5 + tau64*tau64*a4a4 + tau63*tau65*a3a5 + tau65*tau65*a5a5)*12.0;

	//weights
	double s55 = (s1 + s3 + 4.0*s2) / 6.0;
	double s5 = fabs(s6 - s55);
	double r1 = 1.0 + delta*s5 / (s1 + epsilon);
	double r2 = 1.0 + delta*s5 / (s2 + epsilon);
	double r3 = 1.0 + delta*s5 / (s3 + epsilon);
	double r4 = 1.0 + delta*s5 / (s6 + epsilon);
	double a1 = 0.065*r1;
	double a2 = 0.495*r2;
	double a3 = 0.405*r3;
	double a4 = 0.035*r4;
	double tw1 = 1.0 / (a1 + a2 + a3 + a4);
	double w1 = a1*tw1;
	double w2 = a2*tw1;
	double w3 = a3*tw1;
	double w4 = a4*tw1;

	//return weighted average
	return  (w1*(2.0*v1 - 7.0*v2 + 11.0*v3)
		+ w2*(-v2 + 5.0*v3 + 2.0*v4) + w3*(2.0*v3 + 5.0*v4 - v5)
		+ w4*(11.0*v4 - 7.0*v5 + 2.0*v6)) / 6.0;
}
double weno_cu6_M(double *f, double delta)
{

	//assign value to v1, v2,...
	int k = 1;
	double v1 = *(f + k + 2);
	double v2 = *(f + k + 1);
	double v3 = *(f + k);
	double v4 = *(f + k - 1); 
	double v5 = *(f + k - 2);
	double v6 = *(f + k - 3);

	//smoothness indicator
	double epsilon = 1.e-20;
	double s11 = v1 - 2.0*v2 + v3;
	double s12 = v1 - 4.0*v2 + 3.0*v3;
	double s1 = 13.0*s11*s11 + 3.0*s12*s12;
	double s21 = v2 - 2.0*v3 + v4;
	double s22 = v2 - v4;
	double s2 = 13.0*s21*s21 + 3.0*s22*s22;
	double s31 = v3 - 2.0*v4 + v5;
	double s32 = 3.0*v3 - 4.0*v4 + v5;
	double s3 = 13.0*s31*s31 + 3.0*s32*s32;
	double tau61 = (259.0*v6 - 1895.0*v5 + 6670.0*v4 - 2590.0*v3 - 2785.0*v2 + 341.0*v1) / 5760.0;
	double tau62 = -(v5 - 12.0*v4 + 22.0*v3 - 12.0*v2 + v1) / 16.0;
	double tau63 = -(7.0*v6 - 47.0*v5 + 94.0*v4 - 70.0*v3 + 11.0*v2 + 5.0*v1) / 144.0;
	double tau64 = (v5 - 4.0*v4 + 6.0*v3 - 4.0*v2 + v1) / 24.0;
	double tau65 = -(-v6 + 5.0*v5 - 10.0*v4 + 10.0*v3 - 5.0*v2 + v1) / 120.0;
	double a1a1 = 1.0, a2a2 = 13.0 / 3.0, a1a3 = 0.5, a3a3 = 3129.0 / 80.0, a2a4 = 21.0 / 5.0;
	double a1a5 = 1.0 / 8.0, a4a4 = 87617.0 / 140.0, a3a5 = 14127.0 / 224.0, a5a5 = 252337135.0 / 16128.0;
	double s6 = (tau61*tau61*a1a1 + tau62*tau62*a2a2 + tau61*tau63*a1a3 + tau63*tau63*a3a3 + tau62*tau64*a2a4
		+ tau61*tau65*a1a5 + tau64*tau64*a4a4 + tau63*tau65*a3a5 + tau65*tau65*a5a5)*12.0;

	//weights
	double s55 = (s1 + s3 + 4.0*s2) / 6.0;
	double s5 = fabs(s6 - s55);
	double r1 = 1.0 + delta*s5 / (s1 + epsilon);
	double r2 = 1.0 + delta*s5 / (s2 + epsilon);
	double r3 = 1.0 + delta*s5 / (s3 + epsilon);
	double r4 = 1.0 + delta*s5 / (s6 + epsilon);
	double a1 = 0.065*r1;
	double a2 = 0.495*r2;
	double a3 = 0.405*r3;
	double a4 = 0.035*r4;
	double tw1 = 1.0 / (a1 + a2 + a3 + a4);
	double w1 = a1*tw1;
	double w2 = a2*tw1;
	double w3 = a3*tw1;
	double w4 = a4*tw1;

	//return weighted average
	return  (w1*(2.0*v1 - 7.0*v2 + 11.0*v3)
		+ w2*(-v2 + 5.0*v3 + 2.0*v4) + w3*(2.0*v3 + 5.0*v4 - v5)
		+ w4*(11.0*v4 - 7.0*v5 + 2.0*v6)) / 6.0;
}
//-------------------------------------------------------------------------------------------------
// the 6th order WENO Scheme LES
// reference Paper: Hu and Adams,
//Scale separation for implicit large eddy simulation.
//JCP vol 230 (2011) 7240–7249
//-------------------------------------------------------------------------------------------------
double weno_cu6_M1_P(double *f, double delta)
{
	//assign value to v1, v2,...
	int k = 0;
	double v1 = *(f + k - 2);
	double v2 = *(f + k - 1);
	double v3 = *(f + k);
	double v4 = *(f + k + 1); 
	double v5 = *(f + k + 2);
	double v6 = *(f + k + 3);

	//smoothness indicator
	double epsilon = 1.e-3*delta; //1.0e-10;
	double s11 = v1 - 2.0*v2 + v3;
	double s12 = v1 - 4.0*v2 + 3.0*v3;
	double s1 = 13.0*s11*s11 + 3.0*s12*s12;
	double s21 = v2 - 2.0*v3 + v4;
	double s22 = v2 - v4;
	double s2 = 13.0*s21*s21 + 3.0*s22*s22;
	double s31 = v3 - 2.0*v4 + v5;
	double s32 = 3.0*v3 - 4.0*v4 + v5;
	double s3 = 13.0*s31*s31 + 3.0*s32*s32;
	double tau61 = (259.0*v6 - 1895.0*v5 + 6670.0*v4 - 2590.0*v3 - 2785.0*v2 + 341.0*v1)/5760.0;
	double tau62 = - (v5 - 12.0*v4 + 22.0*v3 - 12.0*v2 + v1)/16.0;
	double tau63 = - (7.0*v6 - 47.0*v5 + 94.0*v4 - 70.0*v3 + 11.0*v2 + 5.0*v1)/144.0;
	double tau64 = (v5 - 4.0*v4 + 6.0*v3 - 4.0*v2 + v1)/24.0;
	double tau65 = - (- v6 + 5.0*v5 - 10.0*v4 + 10.0*v3 - 5.0*v2 + v1)/120.0;
	double a1a1 = 1.0, a2a2 = 13.0/3.0, a1a3 = 0.5, a3a3 = 3129.0/80.0, a2a4 = 21.0/5.0;
	double a1a5 = 1.0/8.0, a4a4 = 87617.0/140.0, a3a5 = 14127.0/224.0, a5a5 = 252337135.0/16128.0;
	double s6 = (tau61*tau61*a1a1 + tau62*tau62*a2a2 + tau61*tau63*a1a3 + tau63*tau63*a3a3 + tau62*tau64*a2a4
				  + tau61*tau65*a1a5 + tau64*tau64*a4a4 + tau63*tau65*a3a5 + tau65*tau65*a5a5)*12.0;

	//weights
	double s55 = (s1 + s3 + 4.0*s2)/6.0;
	double s5 = fabs(s6 - s55);
	double r0 = 1.0e3;
	double r1 = r0 + s5/(s1 + epsilon); 
	double r2 = r0 + s5/(s2 + epsilon);
	double r3 = r0 + s5/(s3 + epsilon);
	double r4 = r0 + s5/(s6 + epsilon);
	double a1 = 0.065*r1*r1*r1*r1;
	double a2 = 0.495*r2*r2*r2*r2;
	double a3 = 0.405*r3*r3*r3*r3;
	double a4 = 0.035*r4*r4*r4*r4;
	double tw1= 1.0 / (a1 + a2 + a3 + a4);
	double w1 = a1*tw1;
	double w2 = a2*tw1;
	double w3 = a3*tw1;
	double w4 = a4*tw1;

	//return weighted average
	return  (w1*(2.0*v1 - 7.0*v2 + 11.0*v3)
		  + w2*(-v2 + 5.0*v3 + 2.0*v4) + w3*(2.0*v3 + 5.0*v4 - v5)
		  + w4*(11.0*v4 - 7.0*v5 + 2.0*v6))/6.0;
}
double weno_cu6_M1_M(double *f, double delta)
{

	//assign value to v1, v2,...
	int k = 1;
	double v1 = *(f + k + 2);
	double v2 = *(f + k + 1);
	double v3 = *(f + k);
	double v4 = *(f + k - 1); 
	double v5 = *(f + k - 2);
	double v6 = *(f + k - 3);

	double epsilon = 1.e-3*delta; //1.0e-10;
	double s11 = v1 - 2.0*v2 + v3;
	double s12 = v1 - 4.0*v2 + 3.0*v3;
	double s1 = 13.0*s11*s11 + 3.0*s12*s12;
	double s21 = v2 - 2.0*v3 + v4;
	double s22 = v2 - v4;
	double s2 = 13.0*s21*s21 + 3.0*s22*s22;
	double s31 = v3 - 2.0*v4 + v5;
	double s32 = 3.0*v3 - 4.0*v4 + v5;
	double s3 = 13.0*s31*s31 + 3.0*s32*s32;
	double tau61 = (259.0*v6 - 1895.0*v5 + 6670.0*v4 - 2590.0*v3 - 2785.0*v2 + 341.0*v1)/5760.0;
	double tau62 = - (v5 - 12.0*v4 + 22.0*v3 - 12.0*v2 + v1)/16.0;
	double tau63 = - (7.0*v6 - 47.0*v5 + 94.0*v4 - 70.0*v3 + 11.0*v2 + 5.0*v1)/144.0;
	double tau64 = (v5 - 4.0*v4 + 6.0*v3 - 4.0*v2 + v1)/24.0;
	double tau65 = - (- v6 + 5.0*v5 - 10.0*v4 + 10.0*v3 - 5.0*v2 + v1)/120.0;
	double a1a1 = 1.0, a2a2 = 13.0/3.0, a1a3 = 0.5, a3a3 = 3129.0/80.0, a2a4 = 21.0/5.0;
	double a1a5 = 1.0/8.0, a4a4 = 87617.0/140.0, a3a5 = 14127.0/224.0, a5a5 = 252337135.0/16128.0;
	double s6 = (tau61*tau61*a1a1 + tau62*tau62*a2a2 + tau61*tau63*a1a3 + tau63*tau63*a3a3 + tau62*tau64*a2a4
				  + tau61*tau65*a1a5 + tau64*tau64*a4a4 + tau63*tau65*a3a5 + tau65*tau65*a5a5)*12.0;

	//weights
	double s55 = (s1 + s3 + 4.0*s2)/6.0;
	double s5 = fabs(s6 - s55);
	double r0 = 1.0e3;
	double r1 = r0 + s5/(s1 + epsilon); 
	double r2 = r0 + s5/(s2 + epsilon);
	double r3 = r0 + s5/(s3 + epsilon);
	double r4 = r0 + s5/(s6 + epsilon);
	double a1 = 0.065*r1*r1*r1*r1;
	double a2 = 0.495*r2*r2*r2*r2;
	double a3 = 0.405*r3*r3*r3*r3;
	double a4 = 0.035*r4*r4*r4*r4;
	double tw1= 1.0 / (a1 + a2 + a3 + a4);
	double w1 = a1*tw1;
	double w2 = a2*tw1;
	double w3 = a3*tw1;
	double w4 = a4*tw1;

	//return weighted average
	return  (w1*(2.0*v1 - 7.0*v2 + 11.0*v3)
		  + w2*(-v2 + 5.0*v3 + 2.0*v4) + w3*(2.0*v3 + 5.0*v4 - v5)
		  + w4*(11.0*v4 - 7.0*v5 + 2.0*v6))/6.0;
}
//-------------------------------------------------------------------------------------------------
//		linear scheme for hybrid method
//-------------------------------------------------------------------------------------------------
double du_optimal5(double *f, double delta)
{
	double v1 = *(f - 2);
	double v2 = *(f - 1);
	double v3 = *f;
	double v4 = *(f + 1); 
	double v5 = *(f + 2);
	double v6 = *(f + 3);
	return (v1 - 5.0*v2 + 10.0*v3 - 10.0*v4 + 5.0*v5 - v6)/250.0;
}
double f2_optimal5(double *f, double delta)
{
	double v1 = *(f - 2);
	double v2 = *(f - 1);
	double v3 = *f;
	double v4 = *(f + 1); 
	double v5 = *(f + 2);
	double v6 = *(f + 3);
	return (v1 - 8.0*v2 + 37.0*v3 + 37.0*v4 - 8.0*v5 + v6)/60.0;
}
//-----------------------------------------------------------------------------------------
//		the 2th WENO Scheme
//-----------------------------------------------------------------------------------------
double weno2_P(double *f, double delta)
{
	//assign value to v1, v2,...
	int k = 0;
	double v4 = *(f + k); 
	double v5 = *(f + k + 1);

	//smoothness indicator
	double epsilon = 1.0e-60;
	double s10 = v4;
	double s11 = (v4 - v5);
	double a1a1 = 13.0/12.0;
	double s0 = s10*s10 + epsilon;
	double s1 = s10*s10 + s11*s11*a1a1 + epsilon;
	double tau2 = s11*s11;

	//weights
	double cnst = 1.0;
	double a0 =		   tau2/s0;
	double a1 = cnst + tau2/s1;
	double tw1= 1.0 / (a0 + a1);
	double w0 = a0*tw1;
	double w1 = a1*tw1;

	//return weighted average
	return  w0*v4
		  + w1*0.5*(v4 + v5);
}
double weno2_M(double *f, double delta)
{
	//assign value to v1, v2,...
	int k = 1;
	double v4 = *(f + k); 
	double v5 = *(f + k - 1);

	//smoothness indicator
	double epsilon = 1.0e-60;
	double s10 = v4;
	double s11 = (v4 - v5);
	double a1a1 = 13.0/12.0;
	double s0 = s10*s10 + epsilon;
	double s1 = s10*s10 + s11*s11*a1a1 + epsilon;
	double tau2 = s11*s11;

	//weights
	double cnst = 1.0;
	double a0 =		   tau2/s0;
	double a1 = cnst + tau2/s1;
	double tw1= 1.0 / (a0 + a1);
	double w0 = a0*tw1;
	double w1 = a1*tw1;

	//return weighted average
	return  w0*v4
		  + w1*0.5*(v4 + v5);
}
//-----------------------------------------------------------------------------------------
//		the 3th WENO Scheme
//-----------------------------------------------------------------------------------------
double weno3_P(double *f, double delta)
{
	//assign value to v1, v2,...
	int k = 0;
	double v3 = *(f + k - 1);
	double v4 = *(f + k); 
	double v5 = *(f + k + 1);

	//smoothness indicator
	double epsilon = 1.0e-60;
	double s11 = (v4 - v5);
	double s21 = (v3 - v4);
	double a1a1 = 1.0/12.0;
	double s1 = s11*s11*a1a1 + epsilon;
	double s2 = s21*s21*a1a1 + epsilon;
	double tau3 = (v3 - 2.0*v4 + v5); tau3 = tau3*tau3/320.0;

	//weights
	double cnst = 1.0;
	double a1 = cnst + tau3/s1; a1 = 2.0*a1/3.0;
	double a2 = cnst + tau3/s2; a2 = a2/3.0;
	double tw1= 1.0 / (a1 + a2);
	double w1 = a1*tw1;
	double w2 = a2*tw1;

	//return weighted average
	return  w1*0.5*(v4 + v5)
		  + w2*0.5*(-v3 + 3.0*v4);
}
double weno3_M(double *f, double delta)
{
	//assign value to v1, v2,...
	int k = 1;
	double v3 = *(f + k + 1);
	double v4 = *(f + k); 
	double v5 = *(f + k - 1);

	//smoothness indicator
	double epsilon = 1.0e-60;
	double s11 = (v4 - v5);
	double s21 = (v3 - v4);
	double a1a1 = 1.0/12.0;
	double s1 = s11*s11*a1a1 + epsilon;
	double s2 = s21*s21*a1a1 + epsilon;
	double tau3 = (v3 - 2.0*v4 + v5); tau3 = tau3*tau3/320.0;

	//weights
	double cnst = 1.0;
	double a1 = cnst + tau3/s1; a1 = 2.0*a1/3.0;
	double a2 = cnst + tau3/s2; a2 = a2/3.0;
	double tw1= 1.0 / (a1 + a2);
	double w1 = a1*tw1;
	double w2 = a2*tw1;

	//return weighted average
	return  w1*0.5*(v4 + v5)
		  + w2*0.5*(-v3 + 3.0*v4);
}
//-----------------------------------------------------------------------------------------
//		the 4th WENO Scheme
//-----------------------------------------------------------------------------------------
double weno4_P(double *f, double delta)
{
	//assign value to v1, v2,...
	int k = 0;
	double v3 = *(f + k - 1);
	double v4 = *(f + k); 
	double v5 = *(f + k + 1);
	double v6 = *(f + k + 2);

	double epsilon = 1.0e-8*delta*delta;
	double s11 = (v4 - v5);
	double s21 = (v3 - v4);
	double s31 = (-3.0*v4 + 4.0*v5 - v6)*0.5;
	double s32 = (v4 - 2.0*v5 + v6)*0.5;
	double a1a1 = 1.0, a2a2 = 13.0/3.0, a1a3 = 0.5, a3a3 = 3129.0/80.0, a2a4 = 21.0/5.0;
	double a1a5 = 1.0/8.0, a4a4 = 87617.0/140.0, a3a5 = 14127.0/224.0, a5a5 = 252337135.0/16128.0;
	double s1 = s11*s11*a1a1 + epsilon;
	double s2 = s21*s21*a1a1 + epsilon;
	double s3 = s31*s31*a1a1 + s32*s32*a2a2 + epsilon;
	double tau43 = (-v3 + 3.0*v4 -3.0*v5 + v6)/6.0;
	double tau4 = tau43*tau43*39.0;

	//weights
	double cnst = 4.0e1;
	double a1 = cnst + tau4/s1; a1 = 2.0*a1/6.0;
	double a2 = cnst + tau4/s2; a2 = a2/6.0;
	double a3 = cnst + tau4/s3; a3 = 3.0*a3/6.0;
	double tw1= 1.0 / (a1 + a2 + a3);
	double w1 = a1*tw1;
	double w2 = a2*tw1;
	double w3 = a3*tw1;

	//return weighted average
	return  w1*0.5*(v4 + v5)
		  + w2*0.5*(-v3 + 3.0*v4)
		  + w3*(2.0*v4 + 5.0*v5 - v6)/6.0;
}
double weno4_M(double *f, double delta)
{
	//assign value to v1, v2,...
	int k = 1;
	double v3 = *(f + k + 1);
	double v4 = *(f + k); 
	double v5 = *(f + k - 1);
	double v6 = *(f + k - 2);

	double epsilon = 1.0e-8*delta*delta;
	double s11 = (v4 - v5);
	double s21 = (v3 - v4);
	double s31 = (-3.0*v4 + 4.0*v5 - v6)*0.5;
	double s32 = (v4 - 2.0*v5 + v6)*0.5;
	double a1a1 = 1.0, a2a2 = 13.0/3.0, a1a3 = 0.5, a3a3 = 3129.0/80.0, a2a4 = 21.0/5.0;
	double a1a5 = 1.0/8.0, a4a4 = 87617.0/140.0, a3a5 = 14127.0/224.0, a5a5 = 252337135.0/16128.0;
	double s1 = s11*s11*a1a1 + epsilon;
	double s2 = s21*s21*a1a1 + epsilon;
	double s3 = s31*s31*a1a1 + s32*s32*a2a2 + epsilon;
	double tau43 = (-v3 + 3.0*v4 -3.0*v5 + v6)/6.0;
	double tau4 = tau43*tau43*39.0;

	//weights
	double cnst = 4.0e1;
	double a1 = cnst + tau4/s1; a1 = 2.0*a1/6.0;
	double a2 = cnst + tau4/s2; a2 = a2/6.0;
	double a3 = cnst + tau4/s3; a3 = 3.0*a3/6.0;
	double tw1= 1.0 / (a1 + a2 + a3);
	double w1 = a1*tw1;
	double w2 = a2*tw1;
	double w3 = a3*tw1;

	//return weighted average
	return  w1*0.5*(v4 + v5)
		  + w2*0.5*(-v3 + 3.0*v4)
		  + w3*(2.0*v4 + 5.0*v5 - v6)/6.0;
}
//-----------------------------------------------------------------------------------------
//		the 5th WENO Scheme
//-----------------------------------------------------------------------------------------
double weno5LP_P(double *f, double delta)
{
	//assign value to v1, v2,...
	int k = 0;
	double v1 = *(f + k - 2);
	double v2 = *(f + k - 1);
	double v3 = *(f + k);
	double v4 = *(f + k + 1);
	double v5 = *(f + k + 2);

	//smoothness indicator
	double s1 = (v1 - 2.0*v2 + v3)*(v1 - 2.0*v2 + v3)/720.0
		+ (v1 - 4.0*v2 + 3.0*v3)*(v1 - 4.0*v2 + 3.0*v3)/48.0;
	double s2 = (v2 - 2.0*v3 + v4)*(v2 - 2.0*v3 + v4)/720.0
		+ (v2 - v4)*(v2 - v4)/48.0;
	double s3 = (v3 - 2.0*v4 + v5)*(v3 - 2.0*v4 + v5)/720.0
		+ (3.0*v3 - 4.0*v4 + v5)*(3.0*v3 - 4.0*v4 + v5)/48.0;
	double ds1 = (8.0*v3 - 9.0*v2 + v1)*(8.0*v3 - 9.0*v2 + v1) / 108.0
		+ (2.0*v3 - 3.0*v2 + v1)*(2.0*v3 - 3.0*v2 + v1) / 405.0;
	double ds2 = (v4 - 2.0*v3 + v2)*(v4 - 2.0*v3 + v2) / 45.0
		+ (v2 - v4)*(v2 - v4) / 12.0;
	double ds3 = (v5 - 3.0*v4 + 2.0*v3)*(v5 - 3.0*v4 + 2.0*v3) / 405.0
		+ (v5 - 9.0*v4 + 8.0*v3)*(v5 - 9.0*v4 + 8.0*v3) / 48.0;
	double ss = (s1 + 4.0*s2 + s3)/ 6.0 + 1.0e-6;
	double dss = (ds1 + 4.0*ds2 + ds3) / 6.0 + 1.0e-6;

	//weights
	double a1 = dss*0.1 / ((1.0e-6 + s1)*(1.0e-6 + s1));
	double a2 = dss*0.6 / ((1.0e-6 + s2)*(1.0e-6 + s2));
	double a3 = dss*0.3 / ((1.0e-6 + s3)*(1.0e-6 + s3));
	double da1 = ss / ((1.0e-6 + ds1)*(1.0e-6 + ds1));
	double da2 = ss / ((1.0e-6 + ds2)*(1.0e-6 + ds2));
	double da3 = ss / ((1.0e-6 + ds3)*(1.0e-6 + ds3));

	double tw1 = 1.0 / (a1 + a2 + a3 + da1 + da2 + da3);

	double w1 = a1*tw1;
	double w2 = a2*tw1;
	double w3 = a3*tw1;
	double dw1 = da1*tw1;
	double dw2 = da2*tw1;
	double dw3 = da3*tw1;

	//return weighted average
	return  w1*(2.0*v1 - 7.0*v2 + 11.0*v3) / 6.0
		+ w2*(-v2 + 5.0*v3 + 2.0*v4) / 6.0
		+ w3*(2.0*v3 + 5.0*v4 - v5) / 6.0
		+ dw1*(-v1 + 3.0*v2 + 16.0*v3) / 18.0
		+ dw2*(-v2 + 8.0*v3 - v4) / 6.0
		+ dw3*(16.0*v3 + 3.0*v4 - v5) / 18.0;
}
double weno5LP_M(double *f, double delta)
{
	//assign value to v1, v2,...
	int k = 1;
	double v1 = *(f + k + 2);
	double v2 = *(f + k + 1);
	double v3 = *(f + k);
	double v4 = *(f + k - 1);
	double v5 = *(f + k - 2);

	//smoothness indicator
	double s1 = (v1 - 2.0*v2 + v3)*(v1 - 2.0*v2 + v3) / 720.0
		+ (v1 - 4.0*v2 + 3.0*v3)*(v1 - 4.0*v2 + 3.0*v3) / 48.0;
	double s2 = (v2 - 2.0*v3 + v4)*(v2 - 2.0*v3 + v4) / 720.0
		+ (v2 - v4)*(v2 - v4) / 48.0;
	double s3 = (v3 - 2.0*v4 + v5)*(v3 - 2.0*v4 + v5) / 720.0
		+ (3.0*v3 - 4.0*v4 + v5)*(3.0*v3 - 4.0*v4 + v5) / 48.0;
	double ds1 = (8.0*v3 - 9.0*v2 + v1)*(8.0*v3 - 9.0*v2 + v1) / 108.0
		+ (2.0*v3 - 3.0*v2 + v1)*(2.0*v3 - 3.0*v2 + v1) / 405.0;
	double ds2 = (v4 - 2.0*v3 + v2)*(v4 - 2.0*v3 + v2) / 45.0
		+ (v2 - v4)*(v2 - v4) / 12.0;
	double ds3 = (v5 - 3.0*v4 + 2.0*v3)*(v5 - 3.0*v4 + 2.0*v3) / 405.0
		+ (v5 - 9.0*v4 + 8.0*v3)*(v5 - 9.0*v4 + 8.0*v3) / 48.0;
	double ss = (s1 + 4.0*s2 + s3) / 6.0 + 1.0e-6;
	double dss = (ds1 + 4.0*ds2 + ds3) / 6.0 + 1.0e-6;

	//weights
	double a1 = dss*0.1 / ((1.0e-6 + s1)*(1.0e-6 + s1));
	double a2 = dss*0.6 / ((1.0e-6 + s2)*(1.0e-6 + s2));
	double a3 = dss*0.3 / ((1.0e-6 + s3)*(1.0e-6 + s3));
	double da1 = ss / ((1.0e-6 + ds1)*(1.0e-6 + ds1));
	double da2 = ss / ((1.0e-6 + ds2)*(1.0e-6 + ds2));
	double da3 = ss / ((1.0e-6 + ds3)*(1.0e-6 + ds3));
	
	double tw1 = 1.0 / (a1 + a2 + a3 + da1 + da2 + da3);
	
	double w1 = a1*tw1;
	double w2 = a2*tw1;
	double w3 = a3*tw1;
	double dw1 = da1*tw1;
	double dw2 = da2*tw1;
	double dw3 = da3*tw1;

	//return weighted average
	return  w1*(2.0*v1 - 7.0*v2 + 11.0*v3) / 6.0
		+ w2*(-v2 + 5.0*v3 + 2.0*v4) / 6.0
		+ w3*(2.0*v3 + 5.0*v4 - v5) / 6.0
		+ dw1*(-v1 + 3.0*v2 + 16.0*v3) / 18.0
		+ dw2*(-v2 + 8.0*v3 - v4) / 6.0
		+ dw3*(16.0*v3 + 3.0*v4 - v5) / 18.0;
}
//-----------------------------------------------------------------------------------------
//		the 6th WENO Scheme
//-----------------------------------------------------------------------------------------
double weno6_P(double *f, double delta)
{
	//assign value to v1, v2,...
	int k = 0;
	double v1 = *(f + k - 2);
	double v2 = *(f + k - 1);
	double v3 = *(f + k); 
	double v4 = *(f + k + 1);
	double v5 = *(f + k + 2);
	double v6 = *(f + k + 3);

	//smoothness indicator
	double lambda = (v6 - 5.0*v5 + 10.0*v4 - 10.0*v3 + 5.0*v2 - v1)*(v6 - 5.0*v5 + 10.0*v4 - 10.0*v3 + 5.0*v2 - v1);
	double lambda2 = lambda*lambda;
	double div = 1.0 / (lambda*lambda2*lambda2 + 48703665444.0*lambda2*lambda2 + 1839463528925998080.0*lambda*lambda2
		+ 137354480519857477632000.0*lambda2 + 28304674423347216384000000.0*lambda + 222531556847250309120000000.0);
	double dediv = (772002.0*v6 + 113148.0*v5 + 5334.0*v4 - 795522.0*v3 - 3654.0*v2 - 91308.0*v1)*lambda2*lambda2
		+ (109180288173600.0*v6 + 197594535428640.0*v5 + 30323222375040.0*v4 
		- 1149199336414080.0*v3 + 44137437372000.0*v2 + 767963853064800.0*v1)*lambda*lambda2
		+ (-224743705767102412800.0*v6 + 1391774886950434560000.0*v5 + 323350869425620377600.0*v4
		- 1139934095961383116800.0*v3 - 157426557324138547200.0*v2 - 193021397323430860800.0*v1)*lambda2
		+ (-63879924521422356480000.0*v6 + 222093634169741967360000.0*v5 + 2152242984313193103360000.0*v4
		- 4919558938239060541440000.0*v3 + 2873067133851686338560000.0*v2 - 263964889574138511360000.0*v1)*lambda
		+ (3708859280787505152000000.0*v6 - 29670874246300041216000000.0*v5 + 137227793389137690624000000.0*v4
		- 85303763458112618496000000.0*v3 - 29670874246300041216000000.0*v2 + 3708859280787505152000000.0*v1);

	//return
	return  dediv*div + v3;
}
double weno6_M(double *f, double delta)
{
	//assign value to v1, v2,...
	int k = 1;
	double v1 = *(f + k + 2);
	double v2 = *(f + k + 1);
	double v3 = *(f + k); 
	double v4 = *(f + k - 1);
	double v5 = *(f + k - 2);
	double v6 = *(f + k - 3);

	//smoothness indicator
	double lambda = (v6 - 5.0*v5 + 10.0*v4 - 10.0*v3 + 5.0*v2 - v1)*(v6 - 5.0*v5 + 10.0*v4 - 10.0*v3 + 5.0*v2 - v1);
	double lambda2 = lambda*lambda;
	double div = 1.0 / (lambda*lambda2*lambda2 + 48703665444.0*lambda2*lambda2 + 1839463528925998080.0*lambda*lambda2
		+ 137354480519857477632000.0*lambda2 + 28304674423347216384000000.0*lambda + 222531556847250309120000000.0);
	double dediv = (772002.0*v6 + 113148.0*v5 + 5334.0*v4 - 795522.0*v3 - 3654.0*v2 - 91308.0*v1)*lambda2*lambda2
		+ (109180288173600.0*v6 + 197594535428640.0*v5 + 30323222375040.0*v4
		- 1149199336414080.0*v3 + 44137437372000.0*v2 + 767963853064800.0*v1)*lambda*lambda2
		+ (-224743705767102412800.0*v6 + 1391774886950434560000.0*v5 + 323350869425620377600.0*v4
		- 1139934095961383116800.0*v3 - 157426557324138547200.0*v2 - 193021397323430860800.0*v1)*lambda2
		+ (-63879924521422356480000.0*v6 + 222093634169741967360000.0*v5 + 2152242984313193103360000.0*v4
		- 4919558938239060541440000.0*v3 + 2873067133851686338560000.0*v2 - 263964889574138511360000.0*v1)*lambda
		+ (3708859280787505152000000.0*v6 - 29670874246300041216000000.0*v5 + 137227793389137690624000000.0*v4
		- 85303763458112618496000000.0*v3 - 29670874246300041216000000.0*v2 + 3708859280787505152000000.0*v1);

	//return
	return  dediv*div + v3;
}
//-----------------------------------------------------------------------------------------
//		the 7th WENO Scheme
//-----------------------------------------------------------------------------------------
double weno7_P(double *f, double delta)
{
	//assign value to v1, v2,...
	int k = 0;
	double v1 = *(f + k - 3);
	double v2 = *(f + k - 2);
	double v3 = *(f + k - 1);
	double v4 = *(f + k); 
	double v5 = *(f + k + 1);
	double v6 = *(f + k + 2);
	double v7 = *(f + k + 3);

	//smoothness indicator
	double epsilon = 1.0e-40;
	double s11 = (v4 - v5);
	double s21 = (v3 - v4);
	double s31 = s11;//(-3.0*v4 + 4.0*v5 - v6)*0.5;
	double s32 = (v4 - 2.0*v5 + v6)*0.5;
	double s41 = s21;//(3.0*v4 - 4.0*v3 + v2)*0.5;
	double s42 = (v4 - 2.0*v3 - v2)*0.5;
	double s51 = s11;//(-43.0*v4 + 69.0*v5 - 33.0*v6 + 7.0*v7)/24.0;
	double s52 = s32;//(2.0*v4 - 5.0*v5 + 4.0*v6 - v7)*0.5;
	double s53 = (-v4 + 3.0*v5 - 3.0*v6 + v7)/6.0;
	double s61 = s21;//(43.0*v4 - 69.0*v3 + 33.0*v2 - 7.0*v1)/24.0;
	double s62 = s42;//(2.0*v4 - 5.0*v3 + 4.0*v2 - v1)*0.5;
	double s63 = (v4 - 3.0*v3 + 3.0*v2 - v1)/6.0;
	double a1a1 = 1.0, a2a2 = 13.0/3.0, a1a3 = 0.5, a3a3 = 3129.0/80.0, a2a4 = 21.0/5.0;
	double a1a5 = 1.0/8.0, a4a4 = 87617.0/140.0, a3a5 = 14127.0/224.0, a5a5 = 252337135.0/16128.0;
	double s1 = s11*s11*a1a1 + epsilon;
	double s2 = s21*s21*a1a1 + epsilon;
	double s3 = s31*s31*a1a1 + s32*s32*a2a2 + epsilon;
	double s4 = s41*s41*a1a1 + s42*s42*a2a2 + epsilon;
	double s5 = s51*s51*a1a1 + s52*s52*a2a2 + s51*s53*a1a3 + s53*s53*a3a3 + epsilon;
	double s6 = s61*s61*a1a1 + s62*s62*a2a2 + s61*s63*a1a3 + s63*s63*a3a3 + epsilon;
	double tau71 = -(-259.0*v7 + 2236.0*v6 - 9455.0*v5 + 9455.0*v3 - 2236.0*v2 + 259.0*v1)/11520.0; 
	double tau72 = (37.0*v7 - 462.0*v6 + 3435.0*v5 - 6020.0*v4 + 3435.0*v3 - 462.0*v2 + 37.0*v1)/3840.0; 
	double tau73 = (- 7.0*v7 + 52.0*v6 - 83.0*v5 + 83.0*v3 - 52.0*v2 + 7.0*v1)/288.0; 
	double tau74 = (-5.0*v1 + 54.0*v2 - 171.0*v3 + 244.0*v4 - 171.0*v5 + 54.0*v6 - 5.0*v7)/576.0; 
	double tau75 = (-v1 + 4.0*v2 - 5.0*v3 + 5.0*v5 - 4.0*v6 + v7)/240.0; 
	double tau76 = (v1 - 6.0*v2 + 15.0*v3 - 20.0*v4 + 15.0*v5 - 6.0*v6 + v7)/720.0; 
	double tau7_l = tau71*tau71*a1a1 + tau72*tau72*a2a2 + tau71*tau73*a1a3 + tau73*tau73*a3a3 + epsilon;
	double tau7_h = 576.0*tau74*tau74 + 15600.0*tau75*tau75 + 1440.0*tau74*tau76 + 563220.0*tau76*tau76;
	double tau7 = tau7_h*tau7_l;

	//weights
	double cnst = 1.0e5;
	double a1 = (cnst + tau7/s1); a1 = 12.0*a1*a1*a1*a1;
	double a2 = (cnst + tau7/s2); a2 =  6.0*a2*a2*a2*a2;
	double a3 = (cnst + tau7/s3); a3 =  9.0*a3*a3*a3*a3;
	double a4 = (cnst + tau7/s4); a4 =  3.0*a4*a4*a4*a4;
	double a5 = (cnst + tau7/s5); a5 =  4.0*a5*a5*a5*a5;
	double a6 = (cnst + tau7/s6); a6 =      a6*a6*a6*a6;
	double tw1= 1.0 / (a1 + a2 + a3 + a4 + a5 + a6);
	double w1 = a1*tw1;
	double w2 = a2*tw1;
	double w3 = a3*tw1;
	double w4 = a4*tw1;
	double w5 = a5*tw1;
	double w6 = a6*tw1;

	//return weighted average
	return  w1*0.5*(v4 + v5)
		  + w2*0.5*(-v3 + 3.0*v4)
		  + w3*(2.0*v4 + 5.0*v5 - v6)/6.0
		  + w4*(11.0*v4 - 7.0*v3 + 2.0*v2)/6.0
		  + w5*(3.0*v4 + 13.0*v5 - 5.0*v6 + v7)/12.0
		  + w6*(-3.0*v1 + 13.0*v2 - 23.0*v3 + 25.0*v4)/12.0;
}
double weno7_M(double *f, double delta)
{
	//assign value to v1, v2,...
	int k = 1;
	double v1 = *(f + k + 3);
	double v2 = *(f + k + 2);
	double v3 = *(f + k + 1);
	double v4 = *(f + k); 
	double v5 = *(f + k - 1);
	double v6 = *(f + k - 2);
	double v7 = *(f + k - 3);

	//smoothness indicator
	double epsilon = 1.0e-40;
	double s11 = (v4 - v5);
	double s21 = (v3 - v4);
	double s31 = s11;//(-3.0*v4 + 4.0*v5 - v6)*0.5;
	double s32 = (v4 - 2.0*v5 + v6)*0.5;
	double s41 = s21;//(3.0*v4 - 4.0*v3 + v2)*0.5;
	double s42 = (v4 - 2.0*v3 - v2)*0.5;
	double s51 = s11;//(-43.0*v4 + 69.0*v5 - 33.0*v6 + 7.0*v7)/24.0;
	double s52 = s32;//(2.0*v4 - 5.0*v5 + 4.0*v6 - v7)*0.5;
	double s53 = (-v4 + 3.0*v5 - 3.0*v6 + v7)/6.0;
	double s61 = s21;//(43.0*v4 - 69.0*v3 + 33.0*v2 - 7.0*v1)/24.0;
	double s62 = s42;//(2.0*v4 - 5.0*v3 + 4.0*v2 - v1)*0.5;
	double s63 = (v4 - 3.0*v3 + 3.0*v2 - v1)/6.0;
	double a1a1 = 1.0, a2a2 = 13.0/3.0, a1a3 = 0.5, a3a3 = 3129.0/80.0, a2a4 = 21.0/5.0;
	double a1a5 = 1.0/8.0, a4a4 = 87617.0/140.0, a3a5 = 14127.0/224.0, a5a5 = 252337135.0/16128.0;
	double s1 = s11*s11*a1a1 + epsilon;
	double s2 = s21*s21*a1a1 + epsilon;
	double s3 = s31*s31*a1a1 + s32*s32*a2a2 + epsilon;
	double s4 = s41*s41*a1a1 + s42*s42*a2a2 + epsilon;
	double s5 = s51*s51*a1a1 + s52*s52*a2a2 + s51*s53*a1a3 + s53*s53*a3a3 + epsilon;
	double s6 = s61*s61*a1a1 + s62*s62*a2a2 + s61*s63*a1a3 + s63*s63*a3a3 + epsilon;
	double tau71 = -(-259.0*v7 + 2236.0*v6 - 9455.0*v5 + 9455.0*v3 - 2236.0*v2 + 259.0*v1)/11520.0; 
	double tau72 = (37.0*v7 - 462.0*v6 + 3435.0*v5 - 6020.0*v4 + 3435.0*v3 - 462.0*v2 + 37.0*v1)/3840.0; 
	double tau73 = (- 7.0*v7 + 52.0*v6 - 83.0*v5 + 83.0*v3 - 52.0*v2 + 7.0*v1)/288.0; 
	double tau74 = (-5.0*v1 + 54.0*v2 - 171.0*v3 + 244.0*v4 - 171.0*v5 + 54.0*v6 - 5.0*v7)/576.0; 
	double tau75 = (-v1 + 4.0*v2 - 5.0*v3 + 5.0*v5 - 4.0*v6 + v7)/240.0; 
	double tau76 = (v1 - 6.0*v2 + 15.0*v3 - 20.0*v4 + 15.0*v5 - 6.0*v6 + v7)/720.0; 
	double tau7_l = tau71*tau71*a1a1 + tau72*tau72*a2a2 + tau71*tau73*a1a3 + tau73*tau73*a3a3 + epsilon;
	double tau7_h = 576.0*tau74*tau74 + 15600.0*tau75*tau75 + 1440.0*tau74*tau76 + 563220.0*tau76*tau76;
	double tau7 = tau7_h*tau7_l;

	//weights
	double cnst = 1.0e5;
	double a1 = (cnst + tau7/s1); a1 = 12.0*a1*a1*a1*a1;
	double a2 = (cnst + tau7/s2); a2 =  6.0*a2*a2*a2*a2;
	double a3 = (cnst + tau7/s3); a3 =  9.0*a3*a3*a3*a3;
	double a4 = (cnst + tau7/s4); a4 =  3.0*a4*a4*a4*a4;
	double a5 = (cnst + tau7/s5); a5 =  4.0*a5*a5*a5*a5;
	double a6 = (cnst + tau7/s6); a6 =      a6*a6*a6*a6;
	double tw1= 1.0 / (a1 + a2 + a3 + a4 + a5 + a6);
	double w1 = a1*tw1;
	double w2 = a2*tw1;
	double w3 = a3*tw1;
	double w4 = a4*tw1;
	double w5 = a5*tw1;
	double w6 = a6*tw1;

	//return weighted average
	return  w1*0.5*(v4 + v5)
		  + w2*0.5*(-v3 + 3.0*v4)
		  + w3*(2.0*v4 + 5.0*v5 - v6)/6.0
		  + w4*(11.0*v4 - 7.0*v3 + 2.0*v2)/6.0
		  + w5*(3.0*v4 + 13.0*v5 - 5.0*v6 + v7)/12.0
		  + w6*(-3.0*v1 + 13.0*v2 - 23.0*v3 + 25.0*v4)/12.0;
}
double du_upwind7(double *f, double delta)
{
	//assign value to v1, v2,...
	double v1 = *(f - 3);
	double v2 = *(f - 2);
	double v3 = *(f - 1);
	double v4 = *f; 
	double v5 = *(f + 1);
	double v6 = *(f + 2);
	double v7 = *(f + 3);
	double v8 = *(f + 4);

	return (-3.0*v1 + 21.0*v2 - 63.0*v3 + 105.0*v4 - 105.0*v5 + 63.0*v6 - 21.0*v7 + 3.0*v8)/840.0;
}
double f2_upwind7(double *f, double delta)
{
		//assign value to v1, v2,...
	double v1 = *(f - 3);
	double v2 = *(f - 2);
	double v3 = *(f - 1);
	double v4 = *f; 
	double v5 = *(f + 1);
	double v6 = *(f + 2);
	double v7 = *(f + 3);
	double v8 = *(f + 4);

	return (-3.0*v1 + 29.0*v2 - 139.0*v3 + 533.0*v4 + 533.0*v5 - 139.0*v6 + 29.0*v7 - 3.0*v8)/840.0;
}
//-----------------------------------------------------------------------------------------
//		the 8th WENO Scheme
//-----------------------------------------------------------------------------------------
double weno8_P(double *f, double delta)
{
	//assign value to v1, v2,...
	int k = 0;
	double v1 = *(f + k - 3);
	double v2 = *(f + k - 2);
	double v3 = *(f + k - 1);
	double v4 = *(f + k); 
	double v5 = *(f + k + 1);
	double v6 = *(f + k + 2);
	double v7 = *(f + k + 3);
	double v8 = *(f + k + 4);

	//smoothness indicator
	double epsilon = 1.0e-4*delta*delta;
	double s11 = (v4 - v5);
	double s21 = (v3 - v4);
	double s31 = s11;//(-3.0*v4 + 4.0*v5 - v6)*0.5;
	double s32 = (v4 - 2.0*v5 + v6)*0.5;
	double s41 = s21;//(3.0*v4 - 4.0*v3 + v2)*0.5;
	double s42 = (v4 - 2.0*v3 - v2)*0.5;
	double s51 = s11;//(-43.0*v4 + 69.0*v5 - 33.0*v6 + 7.0*v7)/24.0;
	double s52 = s32;//(2.0*v4 - 5.0*v5 + 4.0*v6 - v7)*0.5;
	double s53 = (-v4 + 3.0*v5 - 3.0*v6 + v7)/6.0;
	double s61 = s21;//(43.0*v4 - 69.0*v3 + 33.0*v2 - 7.0*v1)/24.0;
	double s62 = s42;//(2.0*v4 - 5.0*v3 + 4.0*v2 - v1)*0.5;
	double s63 = (v4 - 3.0*v3 + 3.0*v2 - v1)/6.0;
	double s71 = s11;//(-95.0*v4 + 174.0*v5 - 120.0*v6 + 50.0*v7 - 9.0*v8)/48.0;
	double s72 = s32;//(23.0*v4 - 68.0*v5 + 74.0*v6 - 36.0*v7 + 7.0*v8)/16.0;
	double s73 = s53;//(-5.0*v4 + 18.0*v5 - 24.0*v6 + 14.0*v7 - 3.0*v8)/12.0;
	double s74 = (v4 - 4.0*v5 + 6.0*v6 - 4.0*v7 + v8)/24.0;
	double a1a1 = 1.0, a2a2 = 13.0/3.0, a1a3 = 0.5, a3a3 = 3129.0/80.0, a2a4 = 21.0/5.0;
	double a1a5 = 1.0/8.0, a4a4 = 87617.0/140.0, a3a5 = 14127.0/224.0, a5a5 = 252337135.0/16128.0;
	double s1 = s11*s11*a1a1 + epsilon;
	double s2 = s21*s21*a1a1 + epsilon;
	double s3 = s31*s31*a1a1 + s32*s32*a2a2 + epsilon;
	double s4 = s41*s41*a1a1 + s42*s42*a2a2 + epsilon;
	double s5 = s51*s51*a1a1 + s52*s52*a2a2 + s51*s53*a1a3 + s53*s53*a3a3 + epsilon;
	double s6 = s61*s61*a1a1 + s62*s62*a2a2 + s61*s63*a1a3 + s63*s63*a3a3 + epsilon;
	double s7 = s71*s71*a1a1 + s72*s72*a2a2 + s71*s73*a1a3 + s73*s73*a3a3
			  + s72*s74*a2a4 + s74*s74*a4a4 + epsilon;
	double tau81 = -(3229.0*v8 - 29855.0*v7 + 130417.0*v6 - 377755.0*v5 + 113015.0*v4 + 196931.0*v3 - 40005.0*v2 + 4023.0*v1)/322560.0; 
	double tau82 = (37.0*v7 - 462.0*v6 + 3435.0*v5 - 6020.0*v4 + 3435.0*v3 - 46.0*v2 + 37.0*v1)/3840.0; 
	double tau83 = (141.0*v8 - 1267.0*v7 + 5041.0*v6 - 8255.0*v5 + 4935.0*v4 + 359.0*v3 - 1093.0*v2 + 139.0*v1)/11520.0; 
	double tau84 = -(5.0*v7 - 54.0*v6 + 171.0*v5 - 244.0*v4 + 171.0*v3 - 54.0*v2 + 5.0*v1)/576.0; 
	double tau85 = (-v1 - 5.0*v2 + 43.0*v3 - 105.0*v4 + 125.0*v5 - 79.0*v6 + 25.0*v7 - 3.0*v8)/960.0; 
	double tau86 = (v1 - 6.0*v2 + 15.0*v3 - 20.0*v4 + 15.0*v5 - 6.0*v6 + v7)/720.0; 
	double tau87 = (-v1 + 7.0*v2 - 21.0*v3 + 35.0*v4 - 35.0*v5 + 21.0*v6 - 7.0*v7 + v8)/5040.0; 
	double tau8_l = tau81*tau81 + tau82*tau82*a2a2 + tau81*tau83*a1a3 + tau83*tau83*a3a3
				  + tau82*tau84*a2a4 + tau84*tau84*a4a4 + epsilon;
	double tau8_h = 15600.0*tau85*tau85 + 563220.0*tau86*tau86 + 27599355.0*tau87*tau87;
	double tau8 = tau8_h*tau8_l;

	//weights
	double cnst = 1.0e8;
	double a1 = (cnst + tau8/s1); a1 = 20.0*a1/70.0;
	double a2 = (cnst + tau8/s2); a2 = 10.0*a2/70.0;
	double a3 = (cnst + tau8/s3); a3 = 18.0*a3/70.0;
	double a4 = (cnst + tau8/s4); a4 =  4.0*a4/70.0;
	double a5 = (cnst + tau8/s5); a5 = 12.0*a5/70.0;
	double a6 = (cnst + tau8/s6); a6 =      a6/70.0;
	double a7 = (cnst + tau8/s7); a7 =   5.0*a7/70.0;
	double tw1= 1.0 / (a1 + a2 + a3 + a4 + a5 + a6 + a7);
	double w1 = a1*tw1;
	double w2 = a2*tw1;
	double w3 = a3*tw1;
	double w4 = a4*tw1;
	double w5 = a5*tw1;
	double w6 = a6*tw1;
	double w7 = a7*tw1;

	//return weighted average
	return  w1*0.5*(v4 + v5)
		  + w2*0.5*(-v3 + 3.0*v4)
		  + w3*(2.0*v4 + 5.0*v5 - v6)/6.0
		  + w4*(11.0*v4 - 7.0*v3 + 2.0*v2)/6.0
		  + w5*(3.0*v4 + 13.0*v5 - 5.0*v6 + v7)/12.0
		  + w6*(-3.0*v1 + 13.0*v2 - 23.0*v3 + 25.0*v4)/12.0
		  + w7*(12.0*v4 + 77.0*v5 - 43.0*v6 + 17.0*v7 - 3.0*v8)/60.0;
}
double weno8_M(double *f, double delta)
{
	//assign value to v1, v2,...
	int k = 1;
	double v1 = *(f + k + 3);
	double v2 = *(f + k + 2);
	double v3 = *(f + k + 1);
	double v4 = *(f + k); 
	double v5 = *(f + k - 1);
	double v6 = *(f + k - 2);
	double v7 = *(f + k - 3);
	double v8 = *(f + k - 4);

	//smoothness indicator
	double epsilon = 1.0e-4*delta*delta;
	double s11 = (v4 - v5);
	double s21 = (v3 - v4);
	double s31 = s11;//(-3.0*v4 + 4.0*v5 - v6)*0.5;
	double s32 = (v4 - 2.0*v5 + v6)*0.5;
	double s41 = s21;//(3.0*v4 - 4.0*v3 + v2)*0.5;
	double s42 = (v4 - 2.0*v3 - v2)*0.5;
	double s51 = s11;//(-43.0*v4 + 69.0*v5 - 33.0*v6 + 7.0*v7)/24.0;
	double s52 = s32;//(2.0*v4 - 5.0*v5 + 4.0*v6 - v7)*0.5;
	double s53 = (-v4 + 3.0*v5 - 3.0*v6 + v7)/6.0;
	double s61 = s21;//(43.0*v4 - 69.0*v3 + 33.0*v2 - 7.0*v1)/24.0;
	double s62 = s42;//(2.0*v4 - 5.0*v3 + 4.0*v2 - v1)*0.5;
	double s63 = (v4 - 3.0*v3 + 3.0*v2 - v1)/6.0;
	double s71 = s11;//(-95.0*v4 + 174.0*v5 - 120.0*v6 + 50.0*v7 - 9.0*v8)/48.0;
	double s72 = s32;//(23.0*v4 - 68.0*v5 + 74.0*v6 - 36.0*v7 + 7.0*v8)/16.0;
	double s73 = s53;//(-5.0*v4 + 18.0*v5 - 24.0*v6 + 14.0*v7 - 3.0*v8)/12.0;
	double s74 = (v4 - 4.0*v5 + 6.0*v6 - 4.0*v7 + v8)/24.0;
	double a1a1 = 1.0, a2a2 = 13.0/3.0, a1a3 = 0.5, a3a3 = 3129.0/80.0, a2a4 = 21.0/5.0;
	double a1a5 = 1.0/8.0, a4a4 = 87617.0/140.0, a3a5 = 14127.0/224.0, a5a5 = 252337135.0/16128.0;
	double s1 = s11*s11*a1a1 + epsilon;
	double s2 = s21*s21*a1a1 + epsilon;
	double s3 = s31*s31*a1a1 + s32*s32*a2a2 + epsilon;
	double s4 = s41*s41*a1a1 + s42*s42*a2a2 + epsilon;
	double s5 = s51*s51*a1a1 + s52*s52*a2a2 + s51*s53*a1a3 + s53*s53*a3a3 + epsilon;
	double s6 = s61*s61*a1a1 + s62*s62*a2a2 + s61*s63*a1a3 + s63*s63*a3a3 + epsilon;
	double s7 = s71*s71*a1a1 + s72*s72*a2a2 + s71*s73*a1a3 + s73*s73*a3a3
			  + s72*s74*a2a4 + s74*s74*a4a4 + epsilon;
	double tau81 = -(3229.0*v8 - 29855.0*v7 + 130417.0*v6 - 377755.0*v5 + 113015.0*v4 + 196931.0*v3 - 40005.0*v2 + 4023.0*v1)/322560.0; 
	double tau82 = (37.0*v7 - 462.0*v6 + 3435.0*v5 - 6020.0*v4 + 3435.0*v3 - 46.0*v2 + 37.0*v1)/3840.0; 
	double tau83 = (141.0*v8 - 1267.0*v7 + 5041.0*v6 - 8255.0*v5 + 4935.0*v4 + 359.0*v3 - 1093.0*v2 + 139.0*v1)/11520.0; 
	double tau84 = -(5.0*v7 - 54.0*v6 + 171.0*v5 - 244.0*v4 + 171.0*v3 - 54.0*v2 + 5.0*v1)/576.0; 
	double tau85 = (-v1 - 5.0*v2 + 43.0*v3 - 105.0*v4 + 125.0*v5 - 79.0*v6 + 25.0*v7 - 3.0*v8)/960.0; 
	double tau86 = (v1 - 6.0*v2 + 15.0*v3 - 20.0*v4 + 15.0*v5 - 6.0*v6 + v7)/720.0; 
	double tau87 = (-v1 + 7.0*v2 - 21.0*v3 + 35.0*v4 - 35.0*v5 + 21.0*v6 - 7.0*v7 + v8)/5040.0; 
	double tau8_l = tau81*tau81 + tau82*tau82*a2a2 + tau81*tau83*a1a3 + tau83*tau83*a3a3
				  + tau82*tau84*a2a4 + tau84*tau84*a4a4 + epsilon;
	double tau8_h = 15600.0*tau85*tau85 + 563220.0*tau86*tau86 + 27599355.0*tau87*tau87;
	double tau8 = tau8_h*tau8_l;
	//weights
	double cnst = 1.0e8;
	double a1 = (cnst + tau8/s1); a1 = 20.0*a1/70.0;
	double a2 = (cnst + tau8/s2); a2 = 10.0*a2/70.0;
	double a3 = (cnst + tau8/s3); a3 = 18.0*a3/70.0;
	double a4 = (cnst + tau8/s4); a4 =  4.0*a4/70.0;
	double a5 = (cnst + tau8/s5); a5 = 12.0*a5/70.0;
	double a6 = (cnst + tau8/s6); a6 =      a6/70.0;
	double a7 = (cnst + tau8/s7); a7 =   5.0*a7/70.0;
	double tw1= 1.0 / (a1 + a2 + a3 + a4 + a5 + a6 + a7);
	double w1 = a1*tw1;
	double w2 = a2*tw1;
	double w3 = a3*tw1;
	double w4 = a4*tw1;
	double w5 = a5*tw1;
	double w6 = a6*tw1;
	double w7 = a7*tw1;

	//return weighted average
	return  w1*0.5*(v4 + v5)
		  + w2*0.5*(-v3 + 3.0*v4)
		  + w3*(2.0*v4 + 5.0*v5 - v6)/6.0
		  + w4*(11.0*v4 - 7.0*v3 + 2.0*v2)/6.0
		  + w5*(3.0*v4 + 13.0*v5 - 5.0*v6 + v7)/12.0
		  + w6*(-3.0*v1 + 13.0*v2 - 23.0*v3 + 25.0*v4)/12.0
		  + w7*(12.0*v4 + 77.0*v5 - 43.0*v6 + 17.0*v7 - 3.0*v8)/60.0;
}
//-----------------------------------------------------------------------------------------
//New lower 3rd order WENO 
//-----------------------------------------------------------------------------------------
double weno3old_P(double *f, double delta)
{

	//assign value to v1, v2,...
	int k = 0;
	double v1 = *(f + k - 1);
	double v2 = *(f + k);
	double v3 = *(f + k + 1); 
	
	//smoothness indicator
	double epsilon = 1.0e6*delta*delta;
	double s1 = (v2 - v3)*(v2 - v3) + epsilon; 
	double s2 = (v2 - v1)*(v2 - v1) + epsilon; 

	//weights
	double a1 = 0.666667/s1/s1;
	double a2 = 0.333333/s2/s2;
	double tw1 = 1.0 / (a1 +a2);
	double w1 = tw1*a1;
	double w2 = tw1*a2;

	//return weighted average
	return  w1*(0.5*v2 + 0.5*v3)
		  + w2*(-0.5*v1 + 1.5*v2);
}
double weno3old_M(double *f, double delta)
{

	//assign value to v1, v2,...
	int k = 1;
	double v1 = *(f + k + 1);
	double v2 = *(f + k);
	double v3 = *(f + k - 1); 
	
	//smoothness indicator
	double epsilon = 1.0e6*delta*delta;
	double s1 = (v2 - v3)*(v2 - v3) + epsilon; 
	double s2 = (v2 - v1)*(v2 - v1) + epsilon; 

	//weights
	double a1 = 0.666667/s1/s1;
	double a2 = 0.333333/s2/s2;
	double tw1 = 1.0 / (a1 +a2);
	double w1 = tw1*a1;
	double w2 = tw1*a2;

	//return weighted average
	return  w1*(0.5*v2 + 0.5*v3)
		  + w2*(-0.5*v1 + 1.5*v2);
}
double wenoZ_P(double *f, double delta)
{

	int k;
	double v1, v2, v3, v4, v5;
	double a1, a2, a3, w1, w2, w3;

	//assign value to v1, v2,...
	k = 0;
	v1 = *(f + k - 2);
	v2 = *(f + k - 1);
	v3 = *(f + k);
	v4 = *(f + k + 1); 
	v5 = *(f + k + 2);

	//smoothness indicator
	double epsilon = 1.e-20;
	double s1 = 13.0*(v1 - 2.0*v2 + v3)*(v1 - 2.0*v2 + v3)
		+ 3.0*(v1 - 4.0*v2 + 3.0*v3)*(v1 - 4.0*v2 + 3.0*v3);
	double s2 = 13.0*(v2 - 2.0*v3 + v4)*(v2 - 2.0*v3 + v4)
		+ 3.0*(v2 - v4)*(v2 - v4);
	double s3 = 13.0*(v3 - 2.0*v4 + v5)*(v3 - 2.0*v4 + v5)
		+ 3.0*(3.0*v3 - 4.0*v4 + v5)*(3.0*v3 - 4.0*v4 + v5);
	double ss = fabs(s1 - s3);

	//weights
	a1 = 0.1*(1.0 + delta*ss / (epsilon + s1));
	a2 = 0.6*(1.0 + delta*ss / (epsilon + s2));
	a3 = 0.3*(1.0 + delta*ss / (epsilon + s3));
	double tw1 = 1.0 / (a1 + a2 + a3);
	w1 = a1*tw1;
	w2 = a2*tw1;
	w3 = a3*tw1;

	//return weighted average
	return  w1*(2.0*v1 - 7.0*v2 + 11.0*v3) / 6.0
		+ w2*(-v2 + 5.0*v3 + 2.0*v4) / 6.0
		+ w3*(2.0*v3 + 5.0*v4 - v5) / 6.0;
}
double wenoZ_M(double *f, double delta)
{

	int k;
	double v1, v2, v3, v4, v5;
	double a1, a2,a3, w1, w2, w3;

	//assign value to v1, v2,...
	k = 1;
	v1 = *(f + k + 2);
	v2 = *(f + k + 1);
	v3 = *(f + k);
	v4 = *(f + k - 1); 
	v5 = *(f + k - 2);

	//smoothness indicator
	double epsilon = 1.e-20;
	double s1 = 13.0*(v1 - 2.0*v2 + v3)*(v1 - 2.0*v2 + v3)
		+ 3.0*(v1 - 4.0*v2 + 3.0*v3)*(v1 - 4.0*v2 + 3.0*v3);
	double s2 = 13.0*(v2 - 2.0*v3 + v4)*(v2 - 2.0*v3 + v4)
		+ 3.0*(v2 - v4)*(v2 - v4);
	double s3 = 13.0*(v3 - 2.0*v4 + v5)*(v3 - 2.0*v4 + v5)
		+ 3.0*(3.0*v3 - 4.0*v4 + v5)*(3.0*v3 - 4.0*v4 + v5);
	double ss = fabs(s1 - s3);

	//weights
	a1 = 0.1*(1.0 + delta*ss / (epsilon + s1));
	a2 = 0.6*(1.0 + delta*ss / (epsilon + s2));
	a3 = 0.3*(1.0 + delta*ss / (epsilon + s3));
	double tw1 = 1.0 / (a1 + a2 + a3);
	w1 = a1*tw1;
	w2 = a2*tw1;
	w3 = a3*tw1;

	//return weighted average
	return  w1*(2.0*v1 - 7.0*v2 + 11.0*v3) / 6.0
		+ w2*(-v2 + 5.0*v3 + 2.0*v4) / 6.0
		+ w3*(2.0*v3 + 5.0*v4 - v5) / 6.0;

}
//-----------------------------------------------------------------------------------------
//		fifth order upwind central scheme 
//-----------------------------------------------------------------------------------------
double houc5_P(double *f, double delta)
{
	//assign value to v1, v2,...
	int k = 0;
	double v1 = *(f + k - 2);
	double v2 = *(f + k - 1);
	double v3 = *(f + k);
	double v4 = *(f + k + 1); 
	double v5 = *(f + k + 2);

	//return weighted average
	return  (2.0*v1 -13.0*v2 + 47.0*v3 + 27.0*v4 -3.0*v5)/60.0;

}
double houc5_M(double *f, double delta)
{
	//assign value to v1, v2,...
	int k =	1;
	double v1 = *(f + k + 2);
	double v2 = *(f + k + 1);
	double v3 = *(f + k);
	double v4 = *(f + k - 1); 
	double v5 = *(f + k - 2);

	//return weighted average
	return  (2.0*v1 -13.0*v2 + 47.0*v3 + 27.0*v4 -3.0*v5)/60.0;

}
//-------------------------------------------------------------------------------------------------
//       the 5th WENO Scheme with unifrom accuracy in smooth region
//-------------------------------------------------------------------------------------------------

double weno5u_P(double *f, double delta)
{
	//assign value to v1, v2,...
	int k = 0;
	double v1 = *(f + k - 2);
	double v2 = *(f + k - 1);
	double v3 = *(f + k);
	double v4 = *(f + k + 1); 
	double v5 = *(f + k + 2);
	double v6 = *(f + k + 3);

	//smoothness indicator
	double epsilon = 1.0e-60;
	double s1 = fabs(34.0*v3*v3  - 91.0*v2*v3 + 23.0*v1*v3 + 61.0*v2*v2  - 31.0*v1*v2 + 4.0*v1*v1) + epsilon;
	double s2 = fabs(4.0*v4*v4  - v3*v4 - 7.0*v2*v4 + v3*v3  - v2*v3 + 4.0*v2*v2) + epsilon;
	double s3 = fabs(4.0*v5*v5  - 31.0*v4*v5 + 23.0*v3*v5 + 61.0*v4*v4  - 91.0*v3*v4 + 34.0*v3*v3) + epsilon;
	double s64 = fabs(28847.0*v6*v6  - 427796.0*v5*v6 + 1611952.0*v4*v6 - 576940.0*v3*v6 - 746542.0*v2*v6
			+ 81632.0*v1*v6 + 1589933.0*v5*v5 - 12057644.0*v4*v5 + 4466456.0*v3*v5 + 5436896.0*v2*v5 - 597778.0*v1*v5
			+ 23235668.0*v4*v4  - 18631744.0*v3*v4 - 19576024.0*v2*v4 + 2182124.0*v1*v4 + 5208428.0*v3*v3  
			+ 4953196.0*v2*v3 - 627824.0*v1*v3 + 5546963.0*v2*v2  - 1161452.0*v1*v2 + 61649.0*v1*v1) / 1330560.0 + epsilon; 


	//weights
	double s56 = (11.0*s1 + 11.0*s3 + 38.0*s2)/60.0;
	double s5 = fabs(s64 - s56) + epsilon;
	double delta6 = 1.0;
	double a1 = 0.05*(delta6 + s5/s1);
	double a2 = 0.45*(delta6 + s5/s2);
	double a3 = 0.45*(delta6 + s5/s3);
	double a4 = 0.05*(delta6 + s5/s64);
	double tw1= 1.0 / (a1 + a2 + a3 + a4);
	double w1 = a1*tw1;
	double w2 = a2*tw1;
	double w3 = a3*tw1;
	double w4 = a4*tw1;

	//return weighted average
	return  (w1*(2.0*v1 - 7.0*v2 + 11.0*v3)
		  + w2*(-v2 + 5.0*v3 + 2.0*v4) + w3*(2.0*v3 + 5.0*v4 - v5)
		  + w4*(11.0*v4 - 7.0*v5 + 2.0*v6))/6.0;
}
double weno5u_M(double *f, double delta)
{
	//assign value to v1, v2,...
	int k = 1;
	double v1 = *(f + k + 2);
	double v2 = *(f + k + 1);
	double v3 = *(f + k);
	double v4 = *(f + k - 1); 
	double v5 = *(f + k - 2);
	double v6 = *(f + k - 3);

	//smoothness indicator
	double epsilon = 1.0e-60;
	double s1 = fabs(34.0*v3*v3  - 91.0*v2*v3 + 23.0*v1*v3 + 61.0*v2*v2  - 31.0*v1*v2 + 4.0*v1*v1) + epsilon;
	double s2 = fabs(4.0*v4*v4  - v3*v4 - 7.0*v2*v4 + v3*v3  - v2*v3 + 4.0*v2*v2) + epsilon;
	double s3 = fabs(4.0*v5*v5  - 31.0*v4*v5 + 23.0*v3*v5 + 61.0*v4*v4  - 91.0*v3*v4 + 34.0*v3*v3) + epsilon;
	double s64 = fabs(28847.0*v6*v6  - 427796.0*v5*v6 + 1611952.0*v4*v6 - 576940.0*v3*v6 - 746542.0*v2*v6
			+ 81632.0*v1*v6 + 1589933.0*v5*v5 - 12057644.0*v4*v5 + 4466456.0*v3*v5 + 5436896.0*v2*v5 - 597778.0*v1*v5
			+ 23235668.0*v4*v4  - 18631744.0*v3*v4 - 19576024.0*v2*v4 + 2182124.0*v1*v4 + 5208428.0*v3*v3  
			+ 4953196.0*v2*v3 - 627824.0*v1*v3 + 5546963.0*v2*v2  - 1161452.0*v1*v2 + 61649.0*v1*v1) / 1330560.0 + epsilon; 


	//weights
	double s56 = (11.0*s1 + 11.0*s3 + 38.0*s2)/60.0;
	double s5 = fabs(s64 - s56) + epsilon;
	double delta6 = 1.0;
	double a1 = 0.05*(delta6 + s5/s1);
	double a2 = 0.45*(delta6 + s5/s2);
	double a3 = 0.45*(delta6 + s5/s3);
	double a4 = 0.05*(delta6 + s5/s64);
	double tw1= 1.0 / (a1 + a2 + a3 + a4);
	double w1 = a1*tw1;
	double w2 = a2*tw1;
	double w3 = a3*tw1;
	double w4 = a4*tw1;

	//return weighted average
	return  (w1*(2.0*v1 - 7.0*v2 + 11.0*v3)
		  + w2*(-v2 + 5.0*v3 + 2.0*v4) + w3*(2.0*v3 + 5.0*v4 - v5)
		  + w4*(11.0*v4 - 7.0*v5 + 2.0*v6))/6.0;
}
//-------------------------------------------------------------------------------------------------
//       the 6th WENO Scheme with unifrom accuracy in smooth region
//-------------------------------------------------------------------------------------------------

double weno4u_P(double *f, double delta)
{

	//assign value to v1, v2,...
	int k = 0;
	double v1 = *(f + k - 1);
	double v2 = *(f + k);
	double v3 = *(f + k + 1); 
	double v4 = *(f + k + 2); 
	
	//smoothness indicator
	double epsilon = 1.0e-40;
	double s11 = v2 - v1;
	double s21 = v2 - v3;
	double s31 = v3 - v4;
	double s1 = s11*s11*s11*s11 + epsilon; 
	double s2 = s21*s21*s21*s21 + epsilon; 
	double s3 = s31*s31*s31*s31 + epsilon;
	double s4 = (11.0*s1 - 4.0*s2 + 5.0*s3)/12.0; 
	double s44 = fabs(s1 - s2) + epsilon; 

	//weights
	double cnst = 1.0;
	double a1 = 1.0*(cnst + s44/s1)/6.0;
	double a2 = 4.0*(cnst + s44/s2)/6.0;
	double a3 = 1.0*(cnst + s44/s4)/6.0;
	double w1 = a1/(a1 + a2 + a3);
	double w2 = a2/(a1 + a2 + a3);
	double w3 = a3/(a1 + a2 + a3);


	//return weighted average
	return  w1*(-0.5*v1 + 1.5*v2)
		  + w2*(0.5*v2 + 0.5*v3)
		  + w3*(1.5*v3 - 0.5*v4);
}

double weno4u_M(double *f, double delta)
{

	//assign value to v1, v2,...
	int k = 1;
	double v1 = *(f + k + 1);
	double v2 = *(f + k);
	double v3 = *(f + k - 1); 
	double v4 = *(f + k - 2); 

	//smoothness indicator
	double epsilon = 1.0e-40;
	double s11 = v2 - v1;
	double s21 = v2 - v3;
	double s31 = v3 - v4;
	double s1 = s11*s11*s11*s11 + epsilon; 
	double s2 = s21*s21*s21*s21 + epsilon; 
	double s3 = s31*s31*s31*s31 + epsilon;
	double s4 = (11.0*s1 - 4.0*s2 + 5.0*s3)/12.0; 
	double s44 = fabs(s1 - s2) + epsilon; 

	//weights
	double cnst = 1.0;
	double a1 = 1.0*(cnst + s44/s1)/6.0;
	double a2 = 4.0*(cnst + s44/s2)/6.0;
	double a3 = 1.0*(cnst + s44/s4)/6.0;
	double w1 = a1/(a1 + a2 + a3);
	double w2 = a2/(a1 + a2 + a3);
	double w3 = a3/(a1 + a2 + a3);


	//return weighted average
	return  w1*(-0.5*v1 + 1.5*v2)
		  + w2*(0.5*v2 + 0.5*v3)
		  + w3*(1.5*v3 - 0.5*v4);
}
//-----------------------------------------------------------------------------------------
//		sixth order center scheme
//-----------------------------------------------------------------------------------------
double central_6(double *f, double delta)
{
	//assign value to v1, v2,...
	double v1 = *(f + 3);
	double v2 = *(f + 2);
	double v3 = *(f + 1);
	double v4 = *f; 
	double v5 = *(f - 1);
	double v6 = *(f - 2);

	return (v1 - 8.0*v2 + 37.0*v3 + 37.0*v4 - 8.0*v5 + v6)/60.0;

}//-----------------------------------------------------------------------------------------
//		sixth order center different for derivative
//-----------------------------------------------------------------------------------------
double central_d_6(double *f, double delta)
{
	//assign value to v1, v2,...
	double v1 = *(f + 3);
	double v2 = *(f + 2);
	double v3 = *(f + 1);
	double v4 = *f; 
	double v5 = *(f - 1);
	double v6 = *(f - 2);

	return (2250.0*(v3 - v4) - 125.0*(v2 - v5) + 9.0*(v1 - v6))/1920.0;
//	return (v3 - v4);
//	return (v1 - 9.0*v2 + 45.0*v3 - 45.0*v4 + 9.0*v5 - v6)/60.0;
//	return (- v2 + 9.0*v3 - 9.0*v4 + v5)/12.0;
}
//-----------------------------------------------------------------------------------------
//		resolvable index
//-----------------------------------------------------------------------------------------
double index_r(double *f, double delta)
{
	//assign value to v1, v2,...
	double v4 = *(f + 2);
	double v3 = *(f + 1);
	double v2 = *f; 
	double v1 = *(f - 1);
	double epsilon = 1.0e-40;
	double curv = 0.5*(v3 - 2.0*v2 + v1)*(v3 - 2.0*v2 + v1) + 0.5*(v4 - 2.0*v3 + v2)*(v4 - 2.0*v3 + v2);
	double slope = (v3 - v2)*(v3 - v2);

	return curv /( slope + curv + epsilon);
}
double weno5u_P(double *f, double delta, double r_index)
{
	//assign value to v1, v2,...
	int k = 0;
	double v1 = *(f + k - 2);
	double v2 = *(f + k - 1);
	double v3 = *(f + k);
	double v4 = *(f + k + 1); 
	double v5 = *(f + k + 2);
	double v6 = *(f + k + 3);

	//smoothness indicator
	double epsilon = 1.0e-40;
	double s1 = 13.0*(v1 - 2.0*v2 + v3)*(v1 - 2.0*v2 + v3) 
	   + 3.0*(v1 - 4.0*v2 + 3.0*v3)*(v1 - 4.0*v2 + 3.0*v3) + epsilon;
	double s2 = 13.0*(v2 - 2.0*v3 + v4)*(v2 - 2.0*v3 + v4) 
	   + 3.0*(v2 - v4)*(v2 - v4) + epsilon;
	double s3 = 13.0*(v3 - 2.0*v4 + v5)*(v3 - 2.0*v4 + v5) 
	   + 3.0*(3.0*v3 - 4.0*v4 + v5)*(3.0*v3 - 4.0*v4 + v5) + epsilon;
	double s64 =  fabs(139633.0*v6*v6 - 1429976.0*v5*v6 + 2863984.0*v4*v6 - 2792660.0*v3*v6
			   + 1325006.0*v2*v6 - 245620.0*v1*v6 + 3824847.0*v5*v5  - 15880404.0*v4*v5 + 15929912.0*v3*v5
			   - 7727988.0*v2*v5 + 1458762.0*v1*v5 + 17195652.0*v4*v4  - 35817664.0*v3*v4
			   + 17905032.0*v2*v4 - 3462252.0*v1*v4 + 19510972.0*v3*v3  - 20427884.0*v2*v3
			   + 4086352.0*v1*v3 + 5653317.0*v2*v2  - 2380800.0*v1*v2 + 271779.0*v1*v1) / 10080.0 + epsilon;

	//weights
	double s56 = (s1 + s3 + 4.0*s2)/6.0;
	double s5 = fabs(s64 - s56) + epsilon;
	double delta6 = r_index;
	double a1 = 0.05*(delta6 + s5/s1);
	double a2 = 0.45*(delta6 + s5/s2);
	double a3 = 0.45*(delta6 + s5/s3);
	double a4 = 0.05*(delta6 + s5/s64);
	double tw1= 1.0 / (a1 + a2 + a3 + a4);
	double w1 = a1*tw1;
	double w2 = a2*tw1;
	double w3 = a3*tw1;
	double w4 = a4*tw1;

	//return weighted average
	return  (w1*(2.0*v1 - 7.0*v2 + 11.0*v3)
		  + w2*(-v2 + 5.0*v3 + 2.0*v4) + w3*(2.0*v3 + 5.0*v4 - v5)
		  + w4*(11.0*v4 - 7.0*v5 + 2.0*v6))/6.0;

}
double weno5u_M(double *f, double delta, double r_index)
{

	//assign value to v1, v2,...
	int k = 1;
	double v1 = *(f + k + 2);
	double v2 = *(f + k + 1);
	double v3 = *(f + k);
	double v4 = *(f + k - 1); 
	double v5 = *(f + k - 2);
	double v6 = *(f + k - 3);

	//smoothness indicator
	double epsilon = 1.0e-40;
	double s1 = 13.0*(v1 - 2.0*v2 + v3)*(v1 - 2.0*v2 + v3) 
	   + 3.0*(v1 - 4.0*v2 + 3.0*v3)*(v1 - 4.0*v2 + 3.0*v3) + epsilon;
	double s2 = 13.0*(v2 - 2.0*v3 + v4)*(v2 - 2.0*v3 + v4) 
	   + 3.0*(v2 - v4)*(v2 - v4) + epsilon;
	double s3 = 13.0*(v3 - 2.0*v4 + v5)*(v3 - 2.0*v4 + v5) 
	   + 3.0*(3.0*v3 - 4.0*v4 + v5)*(3.0*v3 - 4.0*v4 + v5) + epsilon;
	double s64 =  fabs(139633.0*v6*v6 - 1429976.0*v5*v6 + 2863984.0*v4*v6 - 2792660.0*v3*v6
			   + 1325006.0*v2*v6 - 245620.0*v1*v6 + 3824847.0*v5*v5  - 15880404.0*v4*v5 + 15929912.0*v3*v5
			   - 7727988.0*v2*v5 + 1458762.0*v1*v5 + 17195652.0*v4*v4  - 35817664.0*v3*v4
			   + 17905032.0*v2*v4 - 3462252.0*v1*v4 + 19510972.0*v3*v3  - 20427884.0*v2*v3
			   + 4086352.0*v1*v3 + 5653317.0*v2*v2  - 2380800.0*v1*v2 + 271779.0*v1*v1) / 10080.0 + epsilon;

	//weights
	double s56 = (s1 + s3 + 4.0*s2)/6.0;
	double s5 = fabs(s64 - s56) + epsilon;
	double delta6 = r_index;
	double a1 = 0.05*(delta6 + s5/s1);
	double a2 = 0.45*(delta6 + s5/s2);
	double a3 = 0.45*(delta6 + s5/s3);
	double a4 = 0.05*(delta6 + s5/s64);
	double tw1= 1.0 / (a1 + a2 + a3 + a4);
	double w1 = a1*tw1;
	double w2 = a2*tw1;
	double w3 = a3*tw1;
	double w4 = a4*tw1;

	//return weighted average
	return  (w1*(2.0*v1 - 7.0*v2 + 11.0*v3)
		  + w2*(-v2 + 5.0*v3 + 2.0*v4) + w3*(2.0*v3 + 5.0*v4 - v5)
		  + w4*(11.0*v4 - 7.0*v5 + 2.0*v6))/6.0;
}
