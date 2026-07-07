// Copyright Avraham Adler (c) 2026
// SPDX-License-Identifier: BSD-2-Clause

#include "Delaporte.h"
#include <R_ext/Rdynload.h>

static const R_CallMethodDef CallEntries[] = {
  {"ddelap_C",    (DL_FUNC) &ddelap_C,   6},
  {"pdelap_C",    (DL_FUNC) &pdelap_C,   7},
  {"qdelap_C",    (DL_FUNC) &qdelap_C,   7},
  {"rdelap_C",    (DL_FUNC) &rdelap_C,   5},
  {"MoMdelap_C",  (DL_FUNC) &MoMdelap_C, 2},
  {"gOMPT_C",     (DL_FUNC) &gOMPT_C,    0},
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
}
