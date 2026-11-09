#ifndef __NBODY_CHK_H__
#define __NBODY_CHK_H__

// 全体系の保存量チェックおよび誤差ダンプ関数
double nbody_chk(const int N, const double eps_th, const double dt,
                 double *XYZUVW, double *KE, double *PE, double *L);

// 従来法（DIM）のバイナリ結果との標準偏差・誤差を比較する関数
void diff_DIM(const double *BIA, const int N);

#endif /* __NBODY_CHK_H__ */