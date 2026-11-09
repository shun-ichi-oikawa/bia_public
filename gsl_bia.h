#ifndef __GSL_BIA_H__
#define __GSL_BIA_H__

#include <gsl/gsl_matrix.h>

// 通常モード用の関数群 (gsl_bia.c 内に実体)
int func(double t, const double *y, double *f, void *params);
int jac(double t, const double *y, double *dfdy, double *dfdt, void *params);

// インクリメンタルモード用の関数群 (__BIA_INC__ 有効時)
int func_inc(double t, const double *Dy, double *f, void *params_inc);
int jac_inc(double t, const double *Dy, double *dfdy, double *dfdt, void *params_inc);

#endif /* __GSL_BIA_H__ */