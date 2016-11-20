#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <string.h>
#include <stdint.h>

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

void qdelap_f(double *p, int np, double *a, int na, double *b, int nb, double *l, int nl,
              int *lt, int *lg, double *ret);


SEXP qdelap_f_wrap(SEXP p, SEXP alpha, SEXP beta, SEXP lambda, SEXP lt, SEXP lg){
  const int np = LENGTH(p);
  const int na = LENGTH(alpha);
  const int nb = LENGTH(beta);
  const int nl = LENGTH(lambda);
  SEXP ret;
  PROTECT(ret = allocVector(REALSXP, np));
  qdelap_f(REAL(p), np, REAL(alpha), na, REAL(beta), nb, REAL(lambda), nl,
           LOGICAL(lt), LOGICAL(lg), REAL(ret));
  UNPROTECT(1);
  return(ret);
}

void rdelap_f(int n, double *a, int na, double *b, int nb, double *l, int nl,
              double *ret);


SEXP rdelap_f_wrap(SEXP n, SEXP alpha, SEXP beta, SEXP lambda){
  const int nn = INTEGER(n)[0];
  const int na = LENGTH(alpha);
  const int nb = LENGTH(beta);
  const int nl = LENGTH(lambda);
  SEXP ret;
  PROTECT(ret = allocVector(REALSXP, nn));
  rdelap_f(nn, REAL(alpha), na, REAL(beta), nb, REAL(lambda), nl, REAL(ret));
  UNPROTECT(1);
  return(ret);
}

void momdelap_f(double *x, int nx, double *ret);

SEXP momdelap_f_wrap(SEXP x){
  const int nx = LENGTH(x);
  SEXP ret;
  PROTECT(ret = allocVector(REALSXP, 3));
  momdelap_f(REAL(x), nx, REAL(ret));
  UNPROTECT(1);
  return(ret);
}

void unifrnd_ (int *n, double *x){
    GetRNGstate();
    for (int i = 0; i < *n; ++i){
        *(x + i) = unif_rand();
    }
    PutRNGstate();
  }

void set_nan_(double *val)
{
  // *val = sqrt(-1.0); By Drew Schmidt
  int64_t x = 0x7FF0000000000001LL;
  memcpy((void *) val, (void *) &x, 8);
}

void set_inf_(double *val) {
  // *val = Inf Based on set_nan
  int64_t x = 0x7FF0000000000000LL;
  memcpy((void *) val, (void *) &x, 8);
}

void set_neginf_(double *val) {
  // *val = Neg Inf Based on set_nan
  int64_t x = 0xFFF0000000000000LL;
  memcpy((void *) val, (void *) &x, 8);
}
