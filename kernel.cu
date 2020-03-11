
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

	return  c/(h*h*h) * tmp1;
	}




int main()
{
	char s[128];
	double dlh, h0, p_wall;
	double x_temp, z_temp;
	int i,j;

	in_file = fopen( "src/Init.txt", "r" );
	fgets(s, 128, in_file);
	fgets(s, 128, in_file);	sscanf(s, "%lf", &Xm);
	fgets(s, 128, in_file);	sscanf(s, "%lf", &Zm);
	fgets(s, 128, in_file);	sscanf(s, "%lf", &X_lq);
	fgets(s, 128, in_file);	sscanf(s, "%lf", &Z_lq);
	fgets(s, 128, in_file);	sscanf(s, "%d", &Ppm);
	fgets(s, 128, in_file); sscanf(s, "%lf", &tau);

	fclose(in_file);

	Pmax = 1000000;

	x = new double[Pmax];
	z = new double[Pmax];
	p = new double[Pmax];
	mas = new double[Pmax];
	rho = new double[Pmax];
	Vx = new double[Pmax];
	Vz = new double[Pmax];
	Ax = new double[Pmax];
	Az = new double[Pmax];
	h = new double[Pmax];
	Ind = new int[Pmax];
	Nn = new int[Pmax];

	dlh = (double) 1.0/Ppm;
	h0 = 4.3*dlh;
	p_wall = 8;

	mas_0 = rho0 / (Ppm);
	mas_0 = mas_0 / (Ppm);

	p_ind = -1;
	Im = int(Xm*Ppm);
	Jm = int(Zm*Ppm);


	// Liquid
	for (i = -Im; i <= Im; i++)
	for (j = 0; j <= Jm; j++)
		{
			x_temp = (double) i*dlh;
			z_temp = (double) j*dlh;
			if ((x_temp <= X_lq) && (z_temp <= Z_lq)){
			p_ind = p_ind + 1;
			x[p_ind] = x_temp;// + (rand()%100-50.0)/100000.0 * dlh;
			z[p_ind] = z_temp;// + (rand()%100-50.0)/100000.0 * dlh;
			rho[p_ind] = rho0;
			p[p_ind] = p0;
			mas[p_ind] = mas_0;
			Vx[p_ind] = 0;
			Vz[p_ind] = 0;
			Az[p_ind] = 0;
			Ind[p_ind] = 0;
			}
		}

	Pr=p_ind; // Last liquid particle

	// Left wall
	for (i = -Im-p_wall; i <= -Im-1; i++)
	for (j = -p_wall; j <= Jm; j++)
		{
			x_temp = (double) i*dlh;
			z_temp = (double) j*dlh;
			p_ind = p_ind + 1;
			x[p_ind] = x_temp;// + (rand()%100-50.0)/100000.0 * dlh;
			z[p_ind] = z_temp;// + (rand()%100-50.0)/100000.0 * dlh;
			rho[p_ind] = rho0;
			p[p_ind] = p0;
			mas[p_ind] = mas_0;
			Vx[p_ind] = 0;
			Vz[p_ind] = 0;
			Az[p_ind] = 0;
			Ind[p_ind] = 1;
		}

	// Right wall
	for (i = Im+1; i <= Im+p_wall; i++)
	for (j = -p_wall; j <= Jm; j++)
		{
			x_temp = (double) i*dlh;
			z_temp = (double) j*dlh;
			p_ind = p_ind + 1;
			x[p_ind] = x_temp;// + (rand()%100-50.0)/100000.0 * dlh;
			z[p_ind] = z_temp;// + (rand()%100-50.0)/100000.0 * dlh;
			rho[p_ind] = rho0;
			p[p_ind] = p0;
			mas[p_ind] = mas_0;
			Vx[p_ind] = 0;
			Vz[p_ind] = 0;
			Az[p_ind] = 0;
			Ind[p_ind] = 1;
		}

	// Bottom
	for (i = -Im; i <= Im; i++)
	for (j = -p_wall; j <= -1; j++)
		{
			x_temp = (double) i*dlh;
			z_temp = (double) j*dlh;
			p_ind = p_ind + 1;
			x[p_ind] = x_temp;// + (rand()%100-50.0)/100000.0 * dlh;
			z[p_ind] = z_temp;// + (rand()%100-50.0)/100000.0 * dlh;
			rho[p_ind] = rho0;
			p[p_ind] = p0;
			mas[p_ind] = mas_0;
			Vx[p_ind] = 0;
			Vz[p_ind] = 0;
			Az[p_ind] = 0;
			Ind[p_ind] = 1;
		}

	Pm=p_ind;




	printf ("Particle: %d \n", Pm);

	out_num = 0;

	out_num = out_num + 1;
	sprintf(out_name, "Data/%d.dat", out_num);
	out_file = fopen( out_name, "w" );
//	fprintf( cut_file, "t=%5.3f mks \n", Tm*1e6 );
	fprintf( out_file, "# x \t z \t rho \t P \t Vx \t Vz \t Ax \t Az \t Ind \t Nn \n" );

    for (p_ind=0; p_ind<=Pm; p_ind++)
    	{
    		fprintf(out_file, "%10.8lf\t%10.8lf\t%10.8lf\t%10.8lf\t%10.8lf\t%10.8lf\t%10.8lf\t%10.8lf\t%d\t%d \n",
    			                          x[p_ind], z[p_ind], rho[p_ind], p[p_ind]/p0, Vx[p_ind], Vz[p_ind], Ax[p_ind], Az[p_ind], Ind[p_ind], Nn[p_ind]);
		}

	fclose(out_file);



	printf ("The end \n");
	return 0;

}
