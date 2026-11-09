//#include <stdio.h>
//#include <math.h>
//#include <gsl/gsl_cblas.h>
void Hermite	( const int, const double, const double  ,       double *,           double *);
void make_aj	( const int, const double, const double *, const double *, double *, double *);

void advance_HMT
	( double *XYZUVW, const double eps_th, const double dt, const int N, const int Nhmt )
{	double r[3*N], v[3*N];
	cblas_dcopy(N,XYZUVW+0,6,r+0,3);cblas_dcopy(N,XYZUVW+1,6,r+1,3);cblas_dcopy(N,XYZUVW+2,6,r+2,3);
	cblas_dcopy(N,XYZUVW+3,6,v+0,3);cblas_dcopy(N,XYZUVW+4,6,v+1,3);cblas_dcopy(N,XYZUVW+5,6,v+2,3);
	for ( int i = 0; i < Nhmt; i++) Hermite( N, eps_th, dt/Nhmt, r, v );
	cblas_dcopy(N,r+0,3,XYZUVW+0,6);cblas_dcopy(N,r+1,3,XYZUVW+1,6);cblas_dcopy(N,r+2,3,XYZUVW+2,6);
	cblas_dcopy(N,v+0,3,XYZUVW+3,6);cblas_dcopy(N,v+1,3,XYZUVW+4,6);cblas_dcopy(N,v+2,3,XYZUVW+5,6);
	return;
}
void Hermite	( const int N, const double eps_th, const double dt, double *r, double *v)
{	double r0[3*N], v0[3*N], a0[3*N], j0[3*N];
	double rp[3*N], vp[3*N], ap[3*N], jp[3*N];
	cblas_dcopy( 3 * N, r, 1, r0, 1 );	/*r0 = r*/
	cblas_dcopy( 3 * N, v, 1, v0, 1 ); /* v0 = v */
	make_aj( N, eps_th, r0, v0, a0, j0 );
	for ( int ip = 0; ip < N; ip++ )	// Predictor
	{	int I = 3 * ip;
		rp[I+0] = r0[I+0] + ( v0[I+0] + ( a0[I+0] + j0[I+0] * dt / 3 ) * dt / 2 ) * dt;
		rp[I+1] = r0[I+1] + ( v0[I+1] + ( a0[I+1] + j0[I+1] * dt / 3 ) * dt / 2 ) * dt;
		rp[I+2] = r0[I+2] + ( v0[I+2] + ( a0[I+2] + j0[I+2] * dt / 3 ) * dt / 2 ) * dt;
		vp[I+0] =             v0[I+0] + ( a0[I+0] + j0[I+0] * dt / 2 ) * dt;
		vp[I+1] =             v0[I+1] + ( a0[I+1] + j0[I+1] * dt / 2 ) * dt;
		vp[I+2] =             v0[I+2] + ( a0[I+2] + j0[I+2] * dt / 2 ) * dt;
	}	make_aj( N, eps_th, rp, vp, ap, jp );
	for ( int ip = 0; ip < N; ip++ )	//	Corrector
	{	int I = 3 * ip;
		v[I+0]  = v0[I+0] + ( a0[I+0] + ap[I+0] + ( j0[I+0] - jp[I+0] ) * dt / 6 ) * dt / 2;
		v[I+1]  = v0[I+1] + ( a0[I+1] + ap[I+1] + ( j0[I+1] - jp[I+1] ) * dt / 6 ) * dt / 2;
		v[I+2]  = v0[I+2] + ( a0[I+2] + ap[I+2] + ( j0[I+2] - jp[I+2] ) * dt / 6 ) * dt / 2;
		r[I+0]  = r0[I+0] + ( v0[I+0] + v [I+0] + ( a0[I+0] - ap[I+0] ) * dt / 6 ) * dt / 2;
		r[I+1]  = r0[I+1] + ( v0[I+1] + v [I+1] + ( a0[I+1] - ap[I+1] ) * dt / 6 ) * dt / 2;
		r[I+2]  = r0[I+2] + ( v0[I+2] + v [I+2] + ( a0[I+2] - ap[I+2] ) * dt / 6 ) * dt / 2;
	}
}
void make_aj
	( const int N, const double eps_th, const double *r, const double *v, double *a, double *j)
{	double r_ik[3], g_ik[3], daik[3], djik[3], dum[3];
	double g_r=0, dist2, dist3, zero[1] = {0};
	cblas_dcopy( 3 * N, zero, 0, a, 1); cblas_dcopy( 3 * N, zero, 0, j, 1);
	//printf("ax  [0]=%11.3e\n",a   [0]);
	for( int ip = 0; ip < N - 1; ip++ )
	{	for(int kp = ip + 1; kp < N; kp++ )
		{	cblas_dcopy( 3,zero,0,daik,1);
			cblas_dcopy( 3,zero,0,djik,1);
			cblas_dcopy( 3,      r + 3 * ip, 1, r_ik, 1 );	//r_ik(:) = r(:,i)
			cblas_daxpy( 3,-1.0, r + 3 * kp, 1, r_ik, 1 );	//        - r(:,k)
			dist2 = cblas_ddot(3,r_ik,1,r_ik,1); dist3 = dist2 * sqrt(dist2);
			//daik(:) = eps * r_ik(:) / dist3;
			cblas_daxpy( 3, eps_th / dist3, r_ik, 1, daik, 1 );
			cblas_daxpy( 3, 1.0, daik, 1, a + 3*ip, 1 );//a(:,i) = a(:,i) + daik;
			cblas_daxpy( 3,-1.0, daik, 1, a + 3*kp, 1 );//a(:,k) = a(:,k) - daik;
			cblas_dcopy( 3,      v + 3 * ip, 1, g_ik, 1 );	//g_ik(:) = v(:,i)
			cblas_daxpy( 3,-1.0, v + 3 * kp, 1, g_ik, 1 );	//        - v(:,k)
			cblas_dcopy( 3,            g_ik, 1,  dum, 1 );
			//djik    = eps / dist3 * ( g_ik(:) - 3 * sum( g_ik * r_ik ) / dist2 * r_ik(:) );
			g_r = cblas_ddot( 3,                       g_ik, 1, r_ik, 1 );
			cblas_daxpy     ( 3, - 3 * g_r    / dist2, r_ik, 1,  dum, 1 );
			cblas_daxpy     ( 3,       eps_th / dist3,  dum, 1, djik, 1 );
			cblas_daxpy( 3, 1.0, djik, 1, j + 3*ip, 1 );	//j(:,i) = j(:,i) + djik(:)
			cblas_daxpy( 3,-1.0, djik, 1, j + 3*kp, 1 );	//j(:,k) = j(:,k) - djik(:)
		}
	}
}
