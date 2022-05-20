/*----------------------------------------------------------------------
  FEMHF_JKLM.h
----------------------------------------------------------------------*/
#ifndef FEMHF_JKLM_H_INCLUDED
#define FEMHF_JKLM_H_INCLUDED

long double FEMHF_K(int k1, int k2, int k3, int k4,
                    int s1, int s2, int s3, int s4, int l);

long double FEMHF_J(int k1, int s1, int k2, int s2, int l);

long double FEMHF_L(int ik, int l, int s1, int s2);

long double FEMHF_M(int ik, int l, int s1, int s2);

void FEMHF_JKLM_Init(int N);
void FEMHF_JKLM_Free(void);

#endif /* FEMHF_JKLM_H_INCLUDED */

