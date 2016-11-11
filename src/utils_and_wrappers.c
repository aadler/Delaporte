#include <R.h>
#include <Rinternals.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

double ddelap_f_s(double *x, double *alpha, double *beta, double *lambda, int *lg);

SEXP ddelap_f_s_wrap(SEXP x, SEXP alpha, SEXP beta, SEXP lambda, SEXP lg){
  SEXP ret;
  PROTECT(ret = allocVector(REALSXP, 1));
  REAL(ret)[0] = ddelap_f_s(REAL(x), REAL(alpha), REAL(beta), REAL(lambda), LOGICAL(lg));
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

SEXP gamln_wrapp(SEXP x){
  SEXP ret;
  PROTECT(ret = allocVector(REALSXP, 1));
  REAL(ret)[0] = gamln(REAL(x));
  UNPROTECT(1);
  return(ret);
}
