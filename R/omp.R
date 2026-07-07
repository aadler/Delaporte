# Copyright (c) 2023, Avraham Adler All rights reserved
# SPDX-License-Identifier: BSD-2-Clause

# Get the number of threads Delaporte will request in its parallel regions.
# This is package-local state (see zzz.R); it does not query---and is not
# affected by---the process-global OpenMP settings of other packages.
getDelapThreads <- function() {
  get("DelapThreads", envir = DelaporteEnv)
}

setDelapThreads <- function(n) {
  n <- as.integer(floor(n[[1L]]))
  if (is.na(n) || n <= 0L) {
    stop("Number of threads must be > 0.")
  }
  ncpus <- get("DLPCPU", envir = DelaporteEnv)
  if (!is.na(ncpus) && n > ncpus) {
    message("Capping at system maximum of ", ncpus, ".")
    n <- ncpus
  }
  assign("DelapThreads", n, envir = DelaporteEnv)
  invisible(n)
}
