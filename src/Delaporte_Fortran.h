// Copyright Avraham Adler (c) 2026
// SPDX-License-Identifier: BSD-2-Clause

#ifndef Delaporte_Fortran_H
#define Delaporte_Fortran_H

extern void gOMPT_f(int *ret);
extern void ddelap_f(double *x, int nx, double *a, int na, double *b, int nb,
                     double *l, int nl, int *lg, int *threads, double *ret);
extern void pdelap_f(double *q, int nq, double *a, int na, double *b, int nb,
                     double *l, int nl, int *lt, int *lg, int *threads,
                     double *ret);

/* First parameter renamed p -> pp to mirror the Fortran dummy, and
   const-qualified as qdelap_f works on an internal copy and never writes to the
  input, so passing REAL(p) directly from the SEXP is now safe.
*/
extern void qdelap_f(const double *pp, int np, double *a, int na, double *b,
                     int nb, double *l, int nl, int *lt, int *lg, int *threads,
                     double *ret);
extern void rdelap_f(int n, double *a, int na, double *b, int nb, double *l,
                     int nl, int *threads, double *ret);
extern void momdelap_f(double *x, int nx, int *tp, double *ret);

#endif
