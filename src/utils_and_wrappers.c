#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <string.h>
#include <stdint.h>

//  Copyright (c) 2016, Avraham Adler
//  All rights reserved.
//  Redistribution and use in source and binary forms, with or without modification, are permitted provided
//  that the following conditions are met:
//    1. Redistributions of source code must retain the above copyright notice, this list of conditions and
//       the following disclaimer.
//    2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
//       and the following disclaimer in the documentation and/or other materials provided with the distribution.
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
//  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
//  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
//  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//  POSSIBILITY OF SUCH DAMAGE.

void ddelap_f(double *x, int nx, double *a, int na, double *b, int nb, double *l, int nl,
              int *lg, double *ret);

SEXP ddelap_C(SEXP x, SEXP alpha, SEXP beta, SEXP lambda, SEXP lg){
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

SEXP pdelap_C(SEXP q, SEXP alpha, SEXP beta, SEXP lambda, SEXP lt, SEXP lg){
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


SEXP qdelap_C(SEXP p, SEXP alpha, SEXP beta, SEXP lambda, SEXP lt, SEXP lg){
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


SEXP rdelap_C(SEXP n, SEXP alpha, SEXP beta, SEXP lambda){
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

SEXP MoMdelap_C(SEXP x){
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
  // *val = sqrt(-1.0); By Drew Schmidt - 2016
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
