#ifndef __IO_H__
#define __IO_H__

#include <stdio.h>

// グローバル変数 EPS の外部参照宣言 (gsl_bia.c等での多重定義エラーを防ぐプロの定石)
extern double EPS;

// 高精度タイマーおよび3次元ベクトル演算
double gettimeofday_sec(void);
double inline dot_product3(double *a, double *b);
void inline cross_product3(const double *a, const double *b, double *c);
void inline cross_product31(const double *a, const double *b, double *e);
double vector3(double *a, double *a_, double *b, double *b_, double *dot, double *crs);

// データ入力・出力インターフェース
void input(double *XYZUVW, const int N);
void input_FORT(double *XYZUVW, const int N);
void output(const double *XYZUVW, const int N);
void out_bia(const double time, const double *XYZUVWY, const int N);
void out_hmt(const double time, const double *XYZUVWY, const int N);
void out_dim(const double time, const double *XYZUVWY, const int N);

// 物理量計算用の内部関数
double energy(const int Np, const double eps_th, const double *XYZUVW,
              double *KE, double *PE, double *LX, double *LY, double *LZ,
              double *UA, double *VA, double *WA);
double Energy(const int Np, const double eps_th, const double *XYZUVW,
              double *KE, double *PE, double *LX, double *LY, double *LZ,
              double *UA, double *VA, double *WA, double *XA, double *YA, double *ZA);

#endif /* __IO_H__ */