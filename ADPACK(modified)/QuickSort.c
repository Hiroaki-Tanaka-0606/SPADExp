/**********************************************************************
  QuickSort.c:

     QuickSort.c is a subroutine to quick-sort an array a with
     an array b. 

  Log of QuickSort.c:

     08/Dec/2005  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adpack.h"

typedef struct {
    double a,b;
} dlists;

typedef struct {
    int a,b;
} ilists;


int dlists_cmp(const dlists *x, const dlists *y);
int ilists_cmp(const ilists *x, const ilists *y);


void qsort_int(long n, int *a, int *b)
{
  int i;
  ilists *AB;

  AB = (ilists*)malloc(sizeof(ilists)*n);

  for (i=0; i<n; i++){
    AB[i].a = a[i+1];     
    AB[i].b = b[i+1];
  }

  qsort(AB, n, sizeof(ilists), (int(*)(const void*, const void*))ilists_cmp);

  for (i=0; i<n; i++){
    a[i+1] = AB[i].a;
    b[i+1] = AB[i].b;
  }

  free(AB);
}


void qsort_double(long n, double *a, double *b)
{
  int i;
  dlists *AB;

  AB = (dlists*)malloc(sizeof(dlists)*n);

  for (i=0; i<n; i++){
    AB[i].a = a[i+1];     
    AB[i].b = b[i+1];
  }

  qsort(AB, n, sizeof(dlists), (int(*)(const void*, const void*))dlists_cmp);

  for (i=0; i<n; i++){
    a[i+1] = AB[i].a;
    b[i+1] = AB[i].b;
  }

  free(AB);
}


 
int dlists_cmp(const dlists *x, const dlists *y)
{
  return (x->a < y->a ? -1 :
          y->a < x->a ?  1 : 0);
}


int ilists_cmp(const ilists *x, const ilists *y)
{
  return (x->a < y->a ? -1 :
          y->a < x->a ?  1 : 0);
}
