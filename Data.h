const double p0=1.0e5;
const double rho0 = 1000.0;
const double Pi = 3.14159265;
const double c0 = 1500.0;

int Pmax;

double Xm, Zm; 				// Full area
double X_lq, Z_lq;			// Liquid area
int Ppm, Pm, Pr, p_ind; 					// Particle per meter
double tau;					// Time step
int Im, Jm;
int out_num;
double mas_0;

double *x,*z;
double *Vx,*Vz;
double *Ax,*Az;
double *mas, *rho, *p;
double *h;
int *Ind, *Nn;


char out_name[25];
FILE *in_file, *out_file, *cut_file;
