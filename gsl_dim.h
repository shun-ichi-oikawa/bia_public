#define METHOD gsl_odeiv2_step_rk8pd
//#define METHOD gsl_odeiv2_step_rkf45
//#define METHOD gsl_odeiv2_step_rk4
//#define METHOD gsl_odeiv2_step_rk2
//#define METHOD gsl_odeiv2_step_rkck
//#define METHOD gsl_odeiv2_step_bsimp
//#define METHOD gsl_odeiv2_step_msbdf
#define GMSET		gsl_matrix_set
const double ABSERR = 1e-14, RELERR = 0;
int func (double t, const double *XYZUVW, double *f, void *params)
{	double eps_th = *(double *)params;
	double zero[1] = {0};
 	double ax[N_P]; cblas_dcopy( N_P, zero, 0, ax, 1);
 	double ay[N_P]; cblas_dcopy( N_P, zero, 0, ay, 1);
 	double az[N_P]; cblas_dcopy( N_P, zero, 0, az, 1);
	for(int i = 0; i < N_P - 1; i++)
	{	int I = 6 * i;
		for(int j = i + 1; j < N_P; j++)
		{	int J = 6 * j;
			double xij = XYZUVW[ I + 0 ] - XYZUVW[ J + 0 ];
			double yij = XYZUVW[ I + 1 ] - XYZUVW[ J + 1 ];
			double zij = XYZUVW[ I + 2 ] - XYZUVW[ J + 2 ];
			double rsq = xij * xij + yij * yij + zij * zij;
			double eee = eps_th / ( rsq * sqrt(rsq) );
			double Fx_ij = xij * eee, Fy_ij = yij * eee, Fz_ij = zij * eee;
			ax[i] += Fx_ij; ax[j] -= Fx_ij;
			ay[i] += Fy_ij; ay[j] -= Fy_ij;
			az[i] += Fz_ij; az[j] -= Fz_ij;
		}
	}
	cblas_dcopy( N_P, XYZUVW + 3, 6, f + 0, 6 ); cblas_dcopy( N_P, ax, 1, f + 3, 6 );
	cblas_dcopy( N_P, XYZUVW + 4, 6, f + 1, 6 ); cblas_dcopy( N_P, ay, 1, f + 4, 6 );
	cblas_dcopy( N_P, XYZUVW + 5, 6, f + 2, 6 ); cblas_dcopy( N_P, az, 1, f + 5, 6 );
	return GSL_SUCCESS;
}

int jac(double t, const double *y, double *dfdy, double *dfdt, void *params)
{	double eps = *(double *)params;
	gsl_matrix_view  dfdy_mat = gsl_matrix_view_array (dfdy, NEQ, NEQ);
	gsl_matrix *m = &dfdy_mat.matrix;
	for(int ip = 0; ip < N_P; ip++)
	{	int IP = 6 * ip;
		for(int jp = 0; jp < N_P; jp++)
		{	int JP = 6 * jp;
			//if( ip == jp )
			//{	for( int i = 0; i < 6; i++)
			//	{	for( int j = 0; j < 6; j++)
			//			GMSET(m,IP+i,JP+j,0.0);
			//	}
			//}
			//else
			{	for(int i = 0; i < 3; i++) /* upper half of 6x6 sub matrix m_ij = dF_i/dY_j */
				{	for(int j = 0; j < 6; j++)
					{	GMSET(m,IP+i,JP+j,(j==i+3)? 1.0 : 0.0); }	//: Z_3x3 + I_3x3 matrix
				}
				for(int i = 3; i < 6; i++) 		/* lower half of 6x6 sub matrix m */
				{	for(int j = 0; j < 3; j++)		/* lower-left of 6x6 sub matrix m */
					{	double r[3] = { y[IP] - y[JP], y[IP+1] - y[JP+1], y[IP+2] - y[JP+2] };
						double r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2], r3 = r2 * sqrt(r2);
						double epr5 = eps / ( r2 * r3 );
						if( ip != jp)
						{	GMSET(m,IP+i,JP+j,(i==j+3)? epr5*(r2-3*r[j]*r[j]): epr5*(-3*r[i-3]*r[j]));
						}
						else
						{	GMSET(m,IP+i,JP+j,0);}
					}
					for(int j = 3; j < 6; j++)		/* lower-right of 6x6 sub matrix m: 3x3 zero matrix */
					{	GMSET(m,IP+i,JP+j,0.0); }
				}
			}
		}
	}
	/* 6x6 matrix of m
	i\j     0 = X        1 = Y        2 = Z        3 = U        4 = V        5 = W
	0 (m,0,0,0.0); (m,0,1,0.0); (m,0,2,0.0); (m,0,3,one); (m,0,4,0.0); (m,0,5,0.0);     _     _
	1 (m,1,0,0.0); (m,1,1,0.0); (m,1,2,0.0); (m,1,3,0.0); (m,1,4,one); (m,1,5,0.0);    | Z   I |
	2 (m,2,0,0.0); (m,2,1,0.0); (m,2,2,0.0); (m,2,3,0.0); (m,2,4,0.0); (m,2,5,one); -> |       |
	3 (m,3,0,m30); (m,3,1,m31); (m,3,2,m32); (m,3,3,0.0); (m,3,4,0.0); (m,3,5,0.0);    |_M   Z_|
	4 (m,4,0,m40); (m,4,1,m41); (m,4,2,m42); (m,4,3,0.0); (m,4,4,0.0); (m,4,5,0.0);
	5 (m,5,0,m50); (m,5,1,m51); (m,5,2,m52); (m,5,3,0.0); (m,5,4,0.0); (m,5,5,0.0);
	lower-left 3x3 sub matrix M of m
	m30=eps*(r2-3*r[3-3]*r[0])/r5,m31=eps*(  -3*r[3-3]*r[1])/r5,m32=eps*(  -3*r[3-3]*r[2])/r5;
	m40=eps*(  -3*r[4-3]*r[0])/r5,m41=eps*(r2-3*r[4-3]*r[1])/r5,m42=eps*(  -3*r[4-3]*r[2])/r5;
	m50=eps*(  -3*r[5-3]*r[0])/r5,m51=eps*(  -3*r[5-3]*r[1])/r5,m52=eps*(r2-3*r[5-3]*r[2])/r5;
	*/
	dfdt[0] = 0.0; dfdt[1] = 0.0; dfdt[2] = 0.0; dfdt[3] = 0.0;
	return GSL_SUCCESS;
}
