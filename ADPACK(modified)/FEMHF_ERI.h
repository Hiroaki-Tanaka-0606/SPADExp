/*----------------------------------------------------------------------
  FEMHF_ERI.h
----------------------------------------------------------------------*/
#ifndef FEMHF_ERI_H_INCLUDED
#define FEMHF_ERI_H_INCLUDED

long double FEMHF_ERI(int k1, int k2, int k3, int k4, 
                      int s1, int s2, int s3, int s4, int l);

long double FEMHF_Gaunt(int l, int m, int l1, int m1, int l2, int m2);

#endif /* FEMHF_ERI_H_INCLUDED */
