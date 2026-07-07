//  Copyright (c) 2016, Avraham Adler
//  All rights reserved.
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided that the following conditions are met:
//    1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions, and the following disclaimer.
//
//    2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
//  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
//  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
//  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
//  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//  POSSIBILITY OF SUCH DAMAGE.



#include "Delaporte.h"
#include "Delaporte_Fortran.h"

/* Reads the OpenMP runtime's max thread count. Used once by .onLoad to seed
 * the package-local thread setting (see R/zzz.R). Its former companion
 * sOMPT_C was removed in 9.0.0: it called omp_set_num_threads(), mutating
 * process-global OpenMP state shared with every other OpenMP-using package.
 * Delaporte now scopes threading per parallel region via num_threads(). */
SEXP gOMPT_C(void) {
  SEXP ret = PROTECT(allocVector(INTSXP, 1));
  gOMPT_f(INTEGER(ret));
  UNPROTECT(1);
  return(ret);
}

SEXP ddelap_C(SEXP x, SEXP alpha, SEXP beta, SEXP lambda, SEXP lg,
              SEXP threads) {
  const int nx = LENGTH(x);
  const int na = LENGTH(alpha);
  const int nb = LENGTH(beta);
  const int nl = LENGTH(lambda);
  /*
   Defense in depth. These entry points are exported via R_RegisterCCallable, so
   external callers can bypass the R-level guards. A zero-length parameter would
   reach imk() in utils.f90 and trigger an integer division by zero (SIGFPE)
   thereby killing the host R process. Now it returns a zero-length real vector,
   mirroring the R wrapper and base R. No PROTECT is needed as the fresh SEXP
   is returned immediately with no intervening allocation.
  */
  if (nx == 0 || na == 0 || nb == 0 || nl == 0) {
    return(allocVector(REALSXP, 0));
  }
  SEXP ret = PROTECT(allocVector(REALSXP, nx));
  ddelap_f(REAL(x), nx, REAL(alpha), na, REAL(beta), nb, REAL(lambda), nl,
           INTEGER(lg), INTEGER(threads), REAL(ret));
  UNPROTECT(1);
  return(ret);
}

SEXP pdelap_C(SEXP q, SEXP alpha, SEXP beta, SEXP lambda, SEXP lt, SEXP lg,
              SEXP threads) {
  const int nq = LENGTH(q);
  const int na = LENGTH(alpha);
  const int nb = LENGTH(beta);
  const int nl = LENGTH(lambda);
  if (nq == 0 || na == 0 || nb == 0 || nl == 0) {
    return(allocVector(REALSXP, 0));
  }
  SEXP ret = PROTECT(allocVector(REALSXP, nq));
  pdelap_f(REAL(q), nq, REAL(alpha), na, REAL(beta), nb, REAL(lambda), nl,
           INTEGER(lt), INTEGER(lg), INTEGER(threads), REAL(ret));
  UNPROTECT(1);
  return(ret);
}

SEXP qdelap_C(SEXP p, SEXP alpha, SEXP beta, SEXP lambda, SEXP lt, SEXP lg,
              SEXP threads) {
  const int np = LENGTH(p);
  const int na = LENGTH(alpha);
  const int nb = LENGTH(beta);
  const int nl = LENGTH(lambda);
  if (np == 0 || na == 0 || nb == 0 || nl == 0) {
    return(allocVector(REALSXP, 0));
  }
  SEXP ret = PROTECT(allocVector(REALSXP, np));
  qdelap_f(REAL(p), np, REAL(alpha), na, REAL(beta), nb, REAL(lambda), nl,
           INTEGER(lt), INTEGER(lg), INTEGER(threads), REAL(ret));
  UNPROTECT(1);
  return(ret);
}

SEXP rdelap_C(SEXP n, SEXP alpha, SEXP beta, SEXP lambda, SEXP threads) {
  const int nn = INTEGER(n)[0];
  const int na = LENGTH(alpha);
  const int nb = LENGTH(beta);
  const int nl = LENGTH(lambda);
  SEXP ret = PROTECT(allocVector(REALSXP, nn));
  /*
   rdelap draws nn variates, so the degenerate-parameter analogue of base R is a
   vector of nn NaNs. The R wrapper attaches the "NaNs produced" warning. Guard
   sits after the allocation so the return object is already the right length.
  */
  if (na == 0 || nb == 0 || nl == 0) {
    double *pret = REAL(ret);
    for (int i = 0; i < nn; ++i) {
      pret[i] = R_NaN;
    }
    UNPROTECT(1);
    return(ret);
  }
  rdelap_f(nn, REAL(alpha), na, REAL(beta), nb, REAL(lambda), nl,
           INTEGER(threads), REAL(ret));
  UNPROTECT(1);
  return(ret);
}

SEXP MoMdelap_C(SEXP x, SEXP tp) {
  const int nx = LENGTH(x);
  SEXP ret = PROTECT(allocVector(REALSXP, 3));
  momdelap_f(REAL(x), nx, INTEGER(tp), REAL(ret));
  UNPROTECT(1);
  return(ret);
}

void unifrnd(const int n, double *x) {
  GetRNGstate();
  for (int i = 0; i < n; ++i) {
    x[i] = unif_rand();
  }
  PutRNGstate();
}
