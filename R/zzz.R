# Copyright (c) 2023, Avraham Adler All rights reserved
# SPDX-License-Identifier: BSD-2-Clause

# Package-local state per Writing R Extensions (r88406): package-specific
# values belong in an environment inside the package namespace---not in the
# global options registry and never in the user's workspace. Defined at
# top level so it is created in the namespace when the package is loaded.

# nolint start: object_name_linter

DelaporteEnv <- new.env(parent = emptyenv())

# nocov start
.onLoad <- function(libname, pkgname) {
  # Get reasonable selection for "max" cpus
  assign("DLPCPU", parallel::detectCores(), envir = DelaporteEnv)
  
  # Initial thread count honors the OpenMP runtime default at load time
  # (e.g. OMP_NUM_THREADS / OMP_THREAD_LIMIT), read once via gOMPT_C. The
  # value is stored package-locally and passed to the Fortran routines on
  # each call; the global OpenMP state is never modified.
  assign("DelapThreads", .Call(gOMPT_C), envir = DelaporteEnv)
}

.onUnload <- function(libpath) {
  # DelaporteEnv lives in the namespace and is garbage-collected with it;
  # nothing to clean up beyond the DLL.
  library.dynam.unload("Delaporte", libpath)
}
# nocov end
# nolint end
