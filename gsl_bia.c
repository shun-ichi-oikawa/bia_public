#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_linalg.h>
//#define __OUTPUT4__
//#define DEBUG
#include "io.h"
double EPS;
#include   "gsl_bia.h"
#include       "bia.h"
#include "nbody_bia.h"
#include "nbody_chk.h"
void advance_HMT( double *, const double, const double, const int );

int main (void)
{	int N; double eps_th;
	FILE *inp; inp = fopen("nbody.inp","r"); fscanf( inp,"%d,%le", &N, &eps_th); fclose(inp);
	double eps = EPS = 2 * eps_th;
	double  *XYZUVW;  XYZUVW = (double*) malloc(6*N*sizeof(double));
	double *dXYZUVW; dXYZUVW = (double*) malloc(6*N*sizeof(double));
	double L_ini[3], L_end[3];	// Need not initialize dXYZUVW = {0}
	//input_FORT(XYZUVW,N);
	input(XYZUVW,N);
	double time = 0, time_end = 1;

	double KE = 0, PE = 0, DE = nbody_chk( N, eps_th, time, XYZUVW, &KE, &PE, L_ini);
	const double E_ini = KE + PE;
	out_bia( time, XYZUVW, 1); //output( XYZUVW, 1);
	double t_start = gettimeofday_sec();
	double t_end_0, t_end_1;
	const int N_time = 1; double dt = time_end / N_time, delta_t = 0;
int Correction = -1;
	for( int i_time = 1; i_time <= N_time; i_time++)
	{	time = i_time * dt; printf("i_time = %d / %d\n", i_time, N_time);
		#ifndef __HMT__
		advance_BIA( XYZUVW, dXYZUVW, eps, dt+delta_t, N); delta_t = 0;
		#else
		advance_HMT( XYZUVW,          eps, dt+delta_t, N); delta_t = 0;
		#endif
		DE = nbody_chk( N, eps_th, time, XYZUVW, &KE, &PE, L_end);
		if	(Correction==1)	//	1: correnction_KE
		{	const double CorFctr = sqrt ( 1 - DE / KE );
			delta_t = ( CorFctr - 1 ) * time;
			correction_KE( N, CorFctr, XYZUVW, delta_t );
			DE = nbody_chk( N, eps_th, time, XYZUVW, &KE, &PE, L_end);
		}
		else if(Correction==0)					//	0: correnction_PE
		{	const double CorFctr =  1 / ( 1 - DE / PE );
			delta_t = ( CorFctr - 1 ) * time;
			correction_PE( N, CorFctr, XYZUVW, delta_t );
			DE = nbody_chk( N, eps_th, time, XYZUVW, &KE, &PE, L_end);
		}
		else if(Correction==-1)
		{	correction_L(N,XYZUVW,L_ini,L_end, E_ini, KE, PE, DE);
			DE = nbody_chk( N, eps_th, time, XYZUVW, &KE, &PE, L_end);
		}
		if ( i_time != N_time ) { out_bia( time, XYZUVW, 1); /*output( XYZUVW, 1);*/ }
	}
	printf("End of time integration.\n\n");
	t_end_0 = gettimeofday_sec();
if(1)
{	diff_DIM( XYZUVW, N);
	advance_BIA( XYZUVW, dXYZUVW, eps, delta_t, N); delta_t = 0;
	DE = nbody_chk( N, eps_th, time, XYZUVW, &KE, &PE, L_end);
	if	( Correction == 1)	//	1: correnction_KE
	{	const double CorFctr = sqrt ( 1 - DE / KE );
		delta_t = ( CorFctr - 1 ) * time;
		correction_KE( N, CorFctr, XYZUVW, delta_t );
	}
	else if( Correction == 0 )					//	0: correnction_PE
	{	const double CorFctr =  1 / ( 1 - DE / PE );
		delta_t = ( CorFctr - 1 ) * time;
		correction_PE( N, CorFctr, XYZUVW, delta_t );
	}
	else
	{	correction_L(N,XYZUVW,L_ini,L_end, E_ini, KE, PE, DE);
		DE = nbody_chk( N, eps_th, time, XYZUVW, &KE, &PE, L_end);
	}
	if ( delta_t != 0)
	{	advance_BIA( XYZUVW, dXYZUVW, eps, delta_t, N); delta_t = 0;
		DE = nbody_chk( N, eps_th, time, XYZUVW, &KE, &PE, L_end);
	}
	diff_DIM( XYZUVW, N);
	for (int itr = 0; itr < 3; itr++)
	{	DE = nbody_chk( N, eps_th, time, XYZUVW, &KE, &PE, L_end);
		correction_L(N,XYZUVW,L_ini,L_end, E_ini, KE, PE, DE);
		DE = nbody_chk( N, eps_th, time, XYZUVW, &KE, &PE, L_end);
		const double lambda = ( sqrt( PE * PE + 4 * KE * E_ini ) - PE ) / ( 2 * KE );
		//printf("lambda = %24.16e, lambda - 1 = %24.16e\n",lambda,lambda-1);
		if(lambda != 1 && Correction == -1 )
		{	cblas_dscal( N, 1 / lambda,XYZUVW+0,6);
			cblas_dscal( N, 1 / lambda,XYZUVW+1,6);
			cblas_dscal( N, 1 / lambda,XYZUVW+2,6);
			cblas_dscal( N,     lambda,XYZUVW+3,6);
			cblas_dscal( N,     lambda,XYZUVW+4,6);
			cblas_dscal( N,     lambda,XYZUVW+5,6);
			correction_L(N,XYZUVW,L_ini,L_end, E_ini, KE, PE, DE);
			DE = nbody_chk( N, eps_th, time, XYZUVW, &KE, &PE, L_end);
		}
	}
}
	out_bia( time, XYZUVW, 1);
	diff_DIM( XYZUVW, N);
	t_end_1 = gettimeofday_sec();
	//printf("CFTR_K=%.16g\n",sqrt(1-DE/KE));
	char filename[20]; sprintf(filename,"%07d.bia",N);
	FILE *bia; bia = fopen( filename,"wb"); fwrite(XYZUVW,sizeof(double),6*N,bia); fclose(bia);
	if(Correction<2)
	{	printf("N = %d, Correction %2d E-time = %.6g sec, %g, %g, DE =%g\n",
			N,Correction,t_end_1-t_end_0,time,time-time_end,DE);
	}
	printf("N = %d, Total         E-time = %.6g sec, %g, %g\n",N,t_end_1-t_start,time,time-time_end);
	free(XYZUVW); free(dXYZUVW);
	return 0;
}
