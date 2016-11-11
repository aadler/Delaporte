#include <R.h>
#include <Rinternals.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

void ddelap_f(double *x, int nx, double *a, int na, double *b, int nb, double *l, int nl, int *lg, double *ret);

SEXP ddelap_f_wrap(SEXP x, SEXP alpha, SEXP beta, SEXP lambda, SEXP lg){
  const int nx = LENGTH(x);
  const int na = LENGTH(alpha);
  const int nb = LENGTH(beta);
  const int nl = LENGTH(lambda);
  SEXP ret;
  PROTECT(ret = allocVector(REALSXP, nx));
  ddelap_f(REAL(x), nx, REAL(alpha), na, REAL(beta), nb, REAL(lambda), nl, LOGICAL(lg), REAL(ret));
  UNPROTECT(1);
  return(ret);
}

void set_nan_(double *val)
{
  // *val = sqrt(-1.0);
  int64_t x = 0x7FF0000000000001LL;
  memcpy((void *) val, (void *) &x, 8);
}

double gamln(double *x);

SEXP gamln_wrap(SEXP x){
  SEXP ret;
  PROTECT(ret = allocVector(REALSXP, 1));
  REAL(ret)[0] = gamln(REAL(x));
  UNPROTECT(1);
  return(ret);
}

void ev_f (double *restrict x, const int nx, double *restrict y, const int ny, double *z);

SEXP ev_f_wrap(SEXP x, SEXP y){
  const int nx = LENGTH(x);
  const int ny = LENGTH(y);
  SEXP z;
  PROTECT(z = allocVector(REALSXP, nx));
  ev_f(REAL(x), nx, REAL(y), ny, REAL(z));
  UNPROTECT(1);
  return(z);
}
