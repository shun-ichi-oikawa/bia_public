void correction_KE
        (       const int N, const double CorFctr, double *XYZUVW, const double delta_t )
{       static int count = 0; count++;
        printf("Begin Correction_K %d\n", count);
                double *UVW;
        for ( int ip = 0; ip < N; ip++)
        {       int I = 6 * ip; UVW = XYZUVW + I + 3;
                {       UVW[0] *= CorFctr; UVW[1] *= CorFctr; UVW[2] *= CorFctr; }
        }
        #ifdef DEBUG
        printf("Current time is changed: dt = %g at correnction_KE\n", delta_t );
        #endif
        printf("End of Correction_K %d, CorFctr=%.16g       \n", count,CorFctr    );
        return;
}

void correction_PE
        (       const int N, const double CorFctr, double *XYZUVW, const double delta_t )
{       static int count = 0; count++;
        printf("Begin Correction_U %d\n", count);
        double *XYZ;
        for ( int ip = 0; ip < N; ip++)
        {       int I = 6 * ip; XYZ = XYZUVW + I;
                XYZ[0] *= CorFctr; XYZ[1] *= CorFctr; XYZ[2] *= CorFctr;
        }
        #ifdef DEBUG
        printf("Current time is changed: dt = %g at correnction_PE\n", delta_t );
        #endif
        printf("End of Correction_U %d, CorFctr=%.16g       \n", count,CorFctr    );
        return;
}

void correction_L       (       const int N, double *XYZUVW, double *L_ini, double *L_end,
                        const double E_ini, const double KE, const double PE, const double DE )
{       static int count = 0; count++;
        printf("Begin Correction_L %d\n", count);

        const double one[1] = {1.0};

        // 1. Save center-of-mass position (R_cm) and velocity (V_cm) in the laboratory frame
        double R_cm[3] = {
                cblas_ddot( N, one, 0, XYZUVW + 0, 6) / N,
                cblas_ddot( N, one, 0, XYZUVW + 1, 6) / N,
                cblas_ddot( N, one, 0, XYZUVW + 2, 6) / N
        };
        double V_cm[3] = {
                cblas_ddot( N, one, 0, XYZUVW + 3, 6) / N,
                cblas_ddot( N, one, 0, XYZUVW + 4, 6) / N,
                cblas_ddot( N, one, 0, XYZUVW + 5, 6) / N
        };

        // 2. Transform to the center-of-mass frame
        for (int d = 0; d < 3; d++) {
                cblas_daxpy( N, -1.0, &R_cm[d], 0, XYZUVW + d,     6);
                cblas_daxpy( N, -1.0, &V_cm[d], 0, XYZUVW + d + 3, 6);
        }

        // Kinetic energy in the center-of-mass frame: KE_cm = KE_lab - 0.5 * N * |V_cm|^2
        const double V_cm_sq = dot_product3(V_cm, V_cm);
        const double KE_cm   = KE - 0.5 * N * V_cm_sq;

        // 3. Coupled angular momentum and energy correction in the center-of-mass frame
        double  A[9] = {0};
        double  r[3], v[3], l[3];
  const double CorFctr = sqrt( 1.0 - DE / KE_cm );
        double dL[3] =  {       L_end[0] - L_ini[0] / CorFctr,
                                                                        L_end[1] - L_ini[1] / CorFctr,
                                                                        L_end[2] - L_ini[2] / CorFctr
                                                                };
        for ( int ip = 0; ip < N; ip++)
        {       int I = 6 * ip;
                r[0] = XYZUVW[I  ]; r[1] = XYZUVW[I+1]; r[2] = XYZUVW[I+2];
                v[0] = XYZUVW[I+3]; v[1] = XYZUVW[I+4]; v[2] = XYZUVW[I+5];
                cross_product3( r, v, l );
                double vsq = dot_product3(v,v);
                if (vsq != 0)
                {       double vsq_inv = 1 / vsq;
                        for (int i=0; i<3; i++)
                        {       for (int j=0; j<3; j++)
                                {       A[3*i+j] += l[i] * v[j] * vsq_inv;
                                }
                        }
                }
        }
        gsl_matrix_view a  = gsl_matrix_view_array(  A, 3, 3);
        gsl_vector_view Dl = gsl_vector_view_array( dL,    3);
        gsl_vector *x = gsl_vector_alloc(3);    // A x = dL will be solved, x = dv
        int s;
        gsl_permutation *p = gsl_permutation_alloc(3);
        gsl_linalg_LU_decomp(&a.matrix, p, &s           );
        gsl_linalg_LU_solve (&a.matrix, p, &Dl.vector, x);
        double dv[3] =
        {       gsl_vector_get ( x, 0 ),
                gsl_vector_get ( x, 1 ),
                gsl_vector_get ( x, 2 ) };
        #ifdef DEBUG
        printf("%5s%12s%12s%12s%15s%15s\n", "","A[][0]", "A[][1]", "A[][2]", "dv[]", "dL[]");
        printf("[0] | %11.6f %11.6f %11.6f  | %11.3e  | %11.3e\n",A[0],A[1],A[2],dv[0],dL[0]);
        printf("[1] | %11.6f %11.6f %11.6f  | %11.3e  | %11.3e\n",A[3],A[4],A[5],dv[1],dL[1]);
        printf("[2] | %11.6f %11.6f %11.6f  | %11.3e  | %11.3e\n",A[6],A[7],A[8],dv[2],dL[2]);
        #endif
        gsl_permutation_free (p);       gsl_vector_free(x);
//        const double CorFctr = sqrt( 1 - DE / KE_cm );
        for ( int ip = 0; ip < N; ip++)
        {       int I = 6 * ip; double *UVW; UVW = XYZUVW + I + 3;
                double V_new[3] = { UVW[0] + dv[0], UVW[1] + dv[1], UVW[2] + dv[2] };
                double V_new_sq = dot_product3(V_new,V_new);
                double v_sq     = dot_product3(  UVW,  UVW);
                for (int i = 0; i < 3; i++) UVW[i] = V_new[i] * sqrt( v_sq / V_new_sq ) * CorFctr;
        }

        // 4. Eliminate residual center-of-mass velocity drift induced by the correction
        double  dU[1] = { cblas_ddot( N, one, 0, XYZUVW + 3, 6) / N };
        double  dV[1] = { cblas_ddot( N, one, 0, XYZUVW + 4, 6) / N };
        double  dW[1] = { cblas_ddot( N, one, 0, XYZUVW + 5, 6) / N };
        cblas_daxpy( N,-1.0, dU, 0, XYZUVW + 3, 6);
        cblas_daxpy( N,-1.0, dV, 0, XYZUVW + 4, 6);
        cblas_daxpy( N,-1.0, dW, 0, XYZUVW + 5, 6);

        // 5. Transform back to the laboratory frame
        for (int d = 0; d < 3; d++) {
                cblas_daxpy( N, +1.0, &R_cm[d], 0, XYZUVW + d,     6);
                cblas_daxpy( N, +1.0, &V_cm[d], 0, XYZUVW + d + 3, 6);
        }

        printf("End of Correction_L %d, CorFctr=%.16g, %.16g\n", count,CorFctr, DE);
        return;
}
