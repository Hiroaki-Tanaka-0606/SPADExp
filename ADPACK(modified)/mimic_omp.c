#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef noomp
#include "mimic_omp.h"

int omp_get_thread_num   (void)
{
  return 0;
}

int omp_get_num_threads  (void)
{
  return 1;
}

int omp_get_num_procs    (void)
{
  return 1;
}

void omp_set_num_threads (int i)
{
}



#endif
