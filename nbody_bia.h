#ifndef __NBODY_BIA_H__
#define __NBODY_BIA_H__

// 運動エネルギー (Kinetic Energy) 補正
void correction_KE(const int N, const double CorFctr, double *XYZUVW, const double delta_t);

// 位置エネルギー (Potential Energy) 補正
void correction_PE(const int N, const double CorFctr, double *XYZUVW, const double delta_t);

// 全角運動量 (Angular Momentum) およびエネルギーの結合補正
void correction_L(const int N, double *XYZUVW, double *L_ini, double *L_end,
                  const double E_ini, const double KE, const double PE, const double DE);

#endif /* __NBODY_BIA_H__ */