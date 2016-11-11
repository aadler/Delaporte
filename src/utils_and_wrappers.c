#include <R.h>
#include <Rinternals.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

void ddelap_f(double *x, int nx, double *a, int na, double *b, int nb, double *l, int nl,
              int *lg, double *ret);

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

void pdelap_f(double *q, int nq, double *a, int na, double *b, int nb, double *l, int nl,
              int *lt, int *lg, double *ret);

SEXP pdelap_f_wrap(SEXP q, SEXP alpha, SEXP beta, SEXP lambda, SEXP lt, SEXP lg){
  const int nq = LENGTH(q);
  const int na = LENGTH(alpha);
  const int nb = LENGTH(beta);
  const int nl = LENGTH(lambda);
  SEXP ret;
  PROTECT(ret = allocVector(REALSXP, nq));
  pdelap_f(REAL(q), nq, REAL(alpha), na, REAL(beta), nb, REAL(lambda), nl,
           LOGICAL(lt), LOGICAL(lg), REAL(ret));
  UNPROTECT(1);
  return(ret);
}

void set_nan_(double *val)
{
  // *val = sqrt(-1.0);
  int64_t x = 0x7FF0000000000001LL;
  memcpy((void *) val, (void *) &x, 8);
}
