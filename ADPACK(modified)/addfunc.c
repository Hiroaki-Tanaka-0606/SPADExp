/**********************************************************************
  addfunc.c:

     addfunc.c is a collective routine of subroutines which are often
     used.

  Log of addfunc.c:

     22/Nov/2001  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adpack.h"

double isgn(int nu)
{
  double result;
  if (nu<0)
    result = -1.0;
  else
    result = 1.0;
  return result;
}

double sgn(double nu)
{
  double result;
  if (nu<0.0)
    result = -1.0;
  else
    result = 1.0;
  return result;
}

double largest(double a, double b)
{
  static double result;

  if (b<=a) result = a;
  else      result = b;
  return result;
}

double smallest(double a, double b)
{
  static double result;

  if (b<=a) result = b;
  else      result = a;
  return result;
}

void fnjoint(char name1[ASIZE8],
             char name2[ASIZE8],
             char name3[ASIZE8])
{
  static char name4[ASIZE8];
  char *f1 = name1,
    *f2 = name2,
    *f3 = name3,
    *f4 = name4;

  while(*f1)
    {
      *f4 = *f1;
      *f1++;
      *f4++;
    }
  while(*f2)
    {
      *f4 = *f2;
      *f2++;
      *f4++;
    }
  while(*f3)
    {
      *f4 = *f3;
      *f3++;
      *f4++;
    }
  *f4 = *f3;
  chcp(name3,name4);
}

void chcp(char name1[ASIZE8],char name2[ASIZE8])
{

  /****************************************************
                    name2 -> name1
  ****************************************************/

  char *f1 = name1,
    *f2 = name2;
  while(*f2){
    *f1 = *f2;
    *f1++;
    *f2++;
  }
  *f1 = *f2;
}
