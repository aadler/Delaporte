// Copyright Avraham Adler (c) 2026
// SPDX-License-Identifier: BSD-2-Clause

#ifndef Delaporte_H
#define Delaporte_H

#include <R.h>
#include <Rinternals.h>

extern SEXP gOMPT_C(void);
extern SEXP ddelap_C(SEXP x, SEXP alpha, SEXP beta, SEXP lambda, SEXP lg,
                     SEXP threads);
extern SEXP pdelap_C(SEXP q, SEXP alpha, SEXP beta, SEXP lambda, SEXP lt,
                     SEXP lg, SEXP threads);
extern SEXP qdelap_C(SEXP p, SEXP alpha, SEXP beta, SEXP lambda, SEXP lt,
                     SEXP lg, SEXP threads);
extern SEXP rdelap_C(SEXP n, SEXP alpha, SEXP beta, SEXP lambda, SEXP threads);
extern SEXP MoMdelap_C(SEXP x, SEXP tp);
void unifrnd(const int n, double *x);

#endif
