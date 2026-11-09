#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
//#define DEBUG
#include "io.h"
#include "nbody_chk.h"
#ifdef __HMT__
#include "gsl_hmt.h"
void advance_HMT	( double *, const double  , const double  , const int   , const int   );
#else
void advance_hmt_	( double *,       double *,       double *,       int * ,       int * );
#endif
int main (void)
{	int N; double eps_th;
	FILE *inp; inp = fopen("nbody.inp","r"); fscanf( inp,"%d,%le", &N, &eps_th); fclose(inp);
	double  *XYZUVW;  XYZUVW = (double*) malloc(6*N*sizeof(double));
	double L_ini[3], L_end[3];
	//input_FORT(XYZUVW,N);
	input(XYZUVW,N);
	double time = 0, time_end = 1;

	double KE = 0, PE = 0, DE = nbody_chk( N, eps_th, time, XYZUVW, &KE, &PE, L_ini);
//	const double E_ini = KE + PE;
	out_hmt( time, XYZUVW, 1); //output( XYZUVW, 1);
	double t_start = gettimeofday_sec();
	double t_end_0, t_end_1;
	const int N_time = 1; double dt = time_end / N_time;
	int Nhmt = 100;
	for( int i_time = 1; i_time <= N_time; i_time++)
	{	time = i_time * dt; printf("i_time = %d / %d\n", i_time, N_time);
		#ifdef __HMT__
		advance_HMT ( XYZUVW,  eps_th,  dt,  N,  Nhmt);
		#else
		advance_hmt_( XYZUVW, &eps_th, &dt, &N, &Nhmt);
		#endif
		DE = nbody_chk( N, eps_th, time, XYZUVW, &KE, &PE, L_end);
		if ( i_time != N_time ) { out_hmt( time, XYZUVW, 1); /*output( XYZUVW, 1);*/ }
	}	//printf("End of time integration.\n");
	t_end_0 = gettimeofday_sec();
	out_hmt( time, XYZUVW, 1);
	diff_DIM( XYZUVW, N);
	t_end_1 = gettimeofday_sec(); char filename[20]; sprintf(filename,"%07d.hmt",N);
	FILE *hmt; hmt = fopen( filename,"wb"); fwrite(XYZUVW,sizeof(double),6*N,hmt); fclose(hmt);
//	printf("N = %d, Correction %2d E-time = %.6g sec, %g, %g, DE =%g\n",
//		N,Correction,t_end_1-t_end_0,time,time-time_end,DE);
	printf("\n");
	printf("N = %d, Total         E-time = %.6g sec, %g, %g\n",N,t_end_1-t_start,time,time-time_end);
	free(XYZUVW);
	return 0;
}
