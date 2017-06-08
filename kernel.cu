
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <cstdlib>
#include <sstream>
#include <string>


#include "Data.h"

//  W, dW  in 3D

/*__device__ __host__ double W(double r, double h)
	{        
	const double Pi = 3.14159265;
	double k, c, tmp1;
	k = (double) fabs(r)/h;
	c= (double) 1/(Pi);
	if (k < 1.0) {tmp1 = (double) 1.0 - 1.5*k*k + 0.75 * k*k*k;}
	if ((k >= 1.0) && (k <= 2.0)) {tmp1 = (double) 0.25 * (2-k)*(2-k)*(2-k);}
	if (k > 2.0) {tmp1 = 0.0;}

	return  c/(h*h*h) * tmp1;
	}

__device__ __host__ double dW(double r, double h)
	{                
	const double Pi = 3.14159265;
	double k, c, tmp1;
	k = (double) r/h;
	c= (double) 1/(Pi);
	if (k < -2.0) {tmp1 = 0.0;}
	if ((k >= -2.0) && (k <= -1.0)) {tmp1 = (double) 0.75 * (2.0+k)*(2.0+k);}
	if ((k > -1.0) && (k < 0)) {tmp1 = (double) -3.0*k - 2.25 * k*k;}
	if ((k >= 0) && (k <= 1.0)) {tmp1 = (double) -3.0*k + 2.25 * k*k;}
	if ((k >= 1.0) && (k <= 2.0)) {tmp1 = (double) -0.75 * (2.0-k)*(2.0-k);}
	if (k > 2.0) {tmp1 = 0.0;}

	return  c/(h*h*h) * tmp1;
	}

*/

//  W,  dW  in 2D


__device__ __host__ double W(double r, double h)
	{
			const double Pi = 3.14159265;
	double k, c, tmp1;
	k = (double) fabs(r)/h;
	c= (double) 10/(7*Pi);
	if (k < 1.0) {tmp1 = (double) 1.0 - 1.5*k*k + 0.75 * k*k*k;}
	if ((k >= 1.0) && (k <= 2.0)) {tmp1 = (double) 0.25 * (2-k)*(2-k)*(2-k);}
	if (k > 2.0) {tmp1 = 0.0;}

	return  c/(h*h) * tmp1;
	}

__device__ __host__ double dW(double r, double h)
	{
			const double Pi = 3.14159265;
	double k, c, tmp1;
	k = (double) r/h;
	c= (double) 10/(7*Pi);
	if (k < -2.0) {tmp1 = 0.0;}
	if ((k >= -2.0) && (k <= -1.0)) {tmp1 = (double) 0.75 * (2.0+k)*(2.0+k);}
	if ((k > -1.0) && (k < 0)) {tmp1 = (double) -3.0*k - 2.25 * k*k;}
	if ((k >= 0) && (k <= 1.0)) {tmp1 = (double) -3.0*k + 2.25 * k*k;}
	if ((k >= 1.0) && (k <= 2.0)) {tmp1 = (double) -0.75 * (2.0-k)*(2.0-k);}
	if (k > 2.0) {tmp1 = 0.0;}

	return  c/(h*h) * tmp1;
	}




int main()
{

	printf ("The end \n");
	return 0;

}
