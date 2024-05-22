#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <Rmath.h>
#include <R_ext/Rdynload.h>

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

void F77_NAME(gOMPT_f)(int *ret);

extern SEXP gOMPT_C(void) {
  SEXP ret = PROTECT(allocVector(INTSXP, 1));
  F77_CALL(gOMPT_f)(INTEGER(ret));
  UNPROTECT(1);
  return(ret);
}

void F77_NAME(sOMPT_f)(int *n);

extern SEXP sOMPT_C(SEXP n) {
  F77_CALL(sOMPT_f)(INTEGER(n));
  return(n);
}

void F77_NAME(ddelap_f)(double *x, int nx, double *a, int na, double *b, 
              int nb, double *l, int nl, int *lg, int *threads, double *ret);

extern SEXP ddelap_C(SEXP x, SEXP alpha, SEXP beta, SEXP lambda, SEXP lg,
                     SEXP threads) {
  const int nx = LENGTH(x);
  const int na = LENGTH(alpha);
  const int nb = LENGTH(beta);
  const int nl = LENGTH(lambda);
  SEXP ret = PROTECT(allocVector(REALSXP, nx));
  F77_CALL(ddelap_f)(REAL(x), nx, REAL(alpha), na, REAL(beta), nb, REAL(lambda),
           nl, INTEGER(lg), INTEGER(threads), REAL(ret));
  UNPROTECT(1);
  return(ret);
}

void F77_NAME(pdelap_f)(double *q, int nq, double *a, int na, double *b, int nb,
              double *l, int nl, int *lt, int *lg, int *threads, double *ret);

extern SEXP pdelap_C(SEXP q, SEXP alpha, SEXP beta, SEXP lambda, SEXP lt,
                     SEXP lg, SEXP threads) {
  const int nq = LENGTH(q);
  const int na = LENGTH(alpha);
  const int nb = LENGTH(beta);
  const int nl = LENGTH(lambda);
  SEXP ret = PROTECT(allocVector(REALSXP, nq));
  F77_CALL(pdelap_f)(REAL(q), nq, REAL(alpha), na, REAL(beta), nb, REAL(lambda),
           nl, INTEGER(lt), INTEGER(lg), INTEGER(threads), REAL(ret));
  UNPROTECT(1);
  return(ret);
}

void F77_NAME(qdelap_f)(double *p, int np, double *a, int na, double *b, int nb,
              double *l, int nl, int *lt, int *lg, int *threads, double *ret);

extern SEXP qdelap_C(SEXP p, SEXP alpha, SEXP beta, SEXP lambda, SEXP lt,
                     SEXP lg, SEXP threads) {
  const int np = LENGTH(p);
  const int na = LENGTH(alpha);
  const int nb = LENGTH(beta);
  const int nl = LENGTH(lambda);
  SEXP ret = PROTECT(allocVector(REALSXP, np));
  F77_CALL(qdelap_f)(REAL(p), np, REAL(alpha), na, REAL(beta), nb, REAL(lambda),
           nl, INTEGER(lt), INTEGER(lg), INTEGER(threads), REAL(ret));
  UNPROTECT(1);
  return(ret);
}

void F77_NAME(rdelap_f)(int n, double *a, int na, double *b, int nb, double *l,
              int nl, int *threads, double *ret);

extern SEXP rdelap_C(SEXP n, SEXP alpha, SEXP beta, SEXP lambda, SEXP threads) {
  const int nn = INTEGER(n)[0];
  const int na = LENGTH(alpha);
  const int nb = LENGTH(beta);
  const int nl = LENGTH(lambda);
  SEXP ret = PROTECT(allocVector(REALSXP, nn));
  F77_CALL(rdelap_f)(nn, REAL(alpha), na, REAL(beta), nb, REAL(lambda), nl,
           INTEGER(threads), REAL(ret));
  UNPROTECT(1);
  return(ret);
}

void F77_NAME(momdelap_f)(double *x, int nx, int *tp, double *ret);

extern SEXP MoMdelap_C(SEXP x, SEXP tp) {
  const int nx = LENGTH(x);
  SEXP ret = PROTECT(allocVector(REALSXP, 3));
  F77_CALL(momdelap_f)(REAL(x), nx, INTEGER(tp), REAL(ret));
  UNPROTECT(1);
  return(ret);
}

void F77_SUB(unifrnd)(int *n, double *x){
  GetRNGstate();
  for (int i = 0; i < *n; ++i){
    *(x + i) = unif_rand();
  }
  PutRNGstate();
}

static const R_CallMethodDef CallEntries[] = {
    {"ddelap_C",    (DL_FUNC) &ddelap_C,   6},
    {"pdelap_C",    (DL_FUNC) &pdelap_C,   7},
    {"qdelap_C",    (DL_FUNC) &qdelap_C,   7},
    {"rdelap_C",    (DL_FUNC) &rdelap_C,   5},
    {"MoMdelap_C",  (DL_FUNC) &MoMdelap_C, 2},
    {"gOMPT_C",     (DL_FUNC) &gOMPT_C,    0},
    {"sOMPT_C",     (DL_FUNC) &sOMPT_C,    1},
    {NULL,                    NULL,        0}
};

void R_init_Delaporte(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
  R_RegisterCCallable("Delaporte", "ddelap_C",  (DL_FUNC) &ddelap_C);
  R_RegisterCCallable("Delaporte", "pdelap_C",  (DL_FUNC) &pdelap_C);
  R_RegisterCCallable("Delaporte", "qdelap_C",  (DL_FUNC) &qdelap_C);
  R_RegisterCCallable("Delaporte", "rdelap_C",  (DL_FUNC) &rdelap_C);
  R_RegisterCCallable("Delaporte", "MoMdelap_C",(DL_FUNC) &MoMdelap_C);
  R_RegisterCCallable("Delaporte", "gOMPT_C",   (DL_FUNC) &gOMPT_C);
  R_RegisterCCallable("Delaporte", "sOMPT_C",   (DL_FUNC) &sOMPT_C);
}
