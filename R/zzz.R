# Copyright (c) 2023, Avraham Adler All rights reserved
# SPDX-License-Identifier: BSD-2-Clause

# nocov start
.onLoad <- function(libname, pkgname) {
  # Get reasonable selection for "max" cpus and store in environment (per 88406)
  DelaporteEnv <<- new.env(parent = emptyenv()) # nolint object_name_linter
  if (!exists("DLPCPU", envir = DelaporteEnv)) {
    assign("DLPCPU", parallel::detectCores(), envir = DelaporteEnv)
  }
}

.onDetach <- function(libpath) {
  # Restore reasonable option
  setDelapThreads(get("DLPCPU", envir = DelaporteEnv))
}

.onUnload <- function(libpath) {
  # Remove holding variable and environment
  rm(DelaporteEnv, pos = 1)
  library.dynam.unload("Delaporte", libpath)
}
# nocov end
