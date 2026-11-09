#include <stdio.h>
//#include <time.h>
//#include <sys/time.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
//#define	DEBUG
#include "io.h"
int N_P;
size_t NEQ;
#include "gsl_dim.h"
#include "nbody_chk.h"
double EPS;

int main (void)
{	double time = 0, time_end = 1;
	double h = ( time_end - time ) / 100;
	FILE *inp; inp = fopen("nbody.inp","r"); fscanf(inp,"%d,%le",&N_P,&EPS); fclose(inp);
	NEQ = 6 * N_P;
	double XYZUVW[NEQ], L_ini[3], L_end[3];
	//input_FORT(XYZUVW,N_P);
	input		(XYZUVW,N_P);

//	gsl_odeiv2_system sys = { func,       NULL  , NEQ, &EPS};
	gsl_odeiv2_system sys = { func, jac /*NULL*/, NEQ, &EPS};
	double KE = 0, PE = 0, DE = nbody_chk( N_P, EPS, time, XYZUVW, &KE, &PE, L_ini);
	//output	(XYZUVW,1); //output(XYZUVW,N_P);
	out_dim	(time,XYZUVW,1);
	double t_start = gettimeofday_sec();
	int itr = 0;
	#ifndef INTERMEDIATE_VALUE
	gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys,METHOD,h,ABSERR,RELERR);
	int status = gsl_odeiv2_driver_apply (d, &time, time_end, XYZUVW);
	if (status != GSL_SUCCESS) {	printf ("GSL error, return value=%d\n", status); }
	out_dim(time,XYZUVW,1);
	gsl_odeiv2_driver_free (d);
	itr++;
	#else
	gsl_odeiv2_step    *s = gsl_odeiv2_step_alloc 		(METHOD, NEQ);
	gsl_odeiv2_control *c = gsl_odeiv2_control_y_new	(ABSERR,RELERR);
	gsl_odeiv2_evolve  *e = gsl_odeiv2_evolve_alloc 	(NEQ);
	while ( time < time_end)
	{	int status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &time, time_end, &h, XYZUVW);
		if (status != GSL_SUCCESS) break;
		out_dim(time,XYZUVW,1); itr++;
	}
	printf ("step method is %s with order ", gsl_odeiv2_step_name (s));
	printf ("%d\n", gsl_odeiv2_step_order(s));
	gsl_odeiv2_evolve_free 	(e);
	gsl_odeiv2_control_free	(c);
	gsl_odeiv2_step_free 	(s);
	#endif
	double t_end = gettimeofday_sec();
	//output(XYZUVW,1);

	char filename[20]; sprintf(filename,"%07d.dim",N_P);
	FILE *out; out = fopen(filename,"wb"); fwrite(XYZUVW,sizeof(double),NEQ,out);fclose(out);

//	correction_L(N_P,XYZUVW,L_ini,L_end);
	DE = nbody_chk( N_P, EPS, time_end, XYZUVW, &KE, &PE, L_end);

	printf("N = %d, iterations = %d to reach time = %g, %g\n", N_P, itr, time, time-time_end);
	printf("Elapsed time = %.6g sec\n", t_end-t_start);
	return 0;
}
