/* shell sort adopted from Edmond Chow */

#include "shellSort_dh.h"

#undef __FUNC__
#define __FUNC__ "shellSort_int"
void shellSort_int(const int n, int *x)
{
  START_FUNC_DH
  int m, max, j, k, itemp;

  m = n/2;
  while (m > 0) {
    max = n - m;
    for (j=0; j<max; j++) {
      for (k=j; k>=0; k-=m) {
        if (x[k+m] >= x[k]) break;
        itemp = x[k+m];
        x[k+m] = x[k];
        x[k] = itemp;
      }
    }
    m = m/2;
  }
  END_FUNC_DH
}


#if 0
#undef __FUNC__
#define __FUNC__ "shellSort_int_float"
void shellSort_int_float(int n, int *x, VAL_DH *xVals)
{
  START_FUNC_DH
  int m, max, j, k, itemp;
  VAL_DH atemp;

  m = n/2;
  while (m > 0) {
    max = n - m;
    for (j=0; j<max; j++) {
      for (k=j; k>=0; k-=m) {
        if (x[k+m] >= x[k]) break;
        itemp = x[k+m];
        atemp = xVals[k+m];
        x[k+m] = x[k];
        /* xVals[k+m] = xVals[k]; */
        x[k] = itemp;
        xVals[k] = atemp;
      }
    }
    m = m/2;
  }
  END_FUNC_DH
}
#endif
