#ifndef __BIA_H__
#define __BIA_H__

// 2体間相対運動のエネルギー評価・出力用関数
double energy4(const double *Y, double *L4);
void output4(const double time, const double *Y);

// BIAアルゴリズムメイン前進サブルーチン
void advance_BIA(double *XYZUVW, double *dXYZUVW, double eps, const double Dt, const int N);

#endif /* __BIA_H__ */