# Copyright (c) 2013, Avraham Adler All rights reserved
# SPDX-License-Identifier: BSD-2-Clause

ddelap <- function(x, alpha, beta, lambda, log = FALSE) {
  # FIX #1: Zero-length parameter vectors previously reached the Fortran
  # recycling helper imk(), where mod(i - 1, 0) is an integer division by zero
  # -> SIGFPE -> the entire R process dies. Guard here and return numeric(0),
  # matching base R semantics (dpois(1, numeric(0)) -> numeric(0)). Zero-length
  # x is also handled explicitly for symmetry and to skip a pointless .Call.
  if (length(alpha) == 0 || length(beta) == 0 || length(lambda) == 0) {
    return(double())
  }
  x <- as.double(x)
  xvalid <- x[!(is.nan(x) | is.na(x))]
  if (any(xvalid > floor(xvalid))) {
    warning("Non-integers passed to ddelap. These will have 0 probability.")
  }
  if (log) log_f <- 1L else log_f <- 0L
  .Call(ddelap_C, x, as.double(alpha), as.double(beta), as.double(lambda),
        log_f, getDelapThreads())
}

pdelap <- function(q, alpha, beta, lambda, lower.tail = TRUE, log.p = FALSE) {
  # Matches R convention to return 0 length result for 0 length parameter. Also
  # Prevents segfault in Fortran.
  if (length(alpha) == 0 || length(beta) == 0 || length(lambda) == 0) {
    return(double())
  }
  
  # These interrupts throw errors even using expect_error. Excluding for now
  # nocov start
  
  # The interactive gate now only guards the genuinely slow routes. With
  # scalar parameters the CDF table is built by an O(K) recurrence
  # (ddelap_table in delaporte.f90), so values up to MAXVECSIZE = 2^24
  # compute in milliseconds and need no gate. What remains quadratic - and
  # therefore still gated at the old 2^15 threshold - is the per-element
  # summation path, taken when any parameter is vector-valued or when q
  # exceeds 2^24. (AA & Claude: 2026-07-07)
  if (any(q[is.finite(q)] >= 2 ^ 63)) {
    stop("Function cannot handle values >= 2^63.")
  }
  
  slowpath <- length(alpha) > 1L || length(beta) > 1L || length(lambda) > 1L
  qmax <- suppressWarnings(max(q[is.finite(q)], -Inf))
  if ((slowpath && qmax >= 2 ^ 15) || qmax >= 2 ^ 24) {
    cat("There are very large values in the supplied data.",
        "This may take minutes if not hours to compute. Are you sure?\n")
    resp <- readline("Press 'y' to continue.\n")
    if (tolower(resp) != "y") {
      cat("Stopping\n")
      return(invisible(NULL))
    }
  }
  # nocov end
  if (lower.tail) lt_f <- 1L else lt_f <- 0L
  if (log.p) lp_f <- 1L else lp_f <- 0L
  .Call(pdelap_C, as.double(q), as.double(alpha), as.double(beta),
        as.double(lambda), lt_f, lp_f, getDelapThreads())
}

qdelap <- function(p, alpha, beta, lambda, lower.tail = TRUE, log.p = FALSE,
                   exact = TRUE) {
  # Matches R convention to return 0 length result for 0 length parameter. Also
  # Prevents segfault in Fortran.
  if (length(alpha) == 0 || length(beta) == 0 || length(lambda) == 0) {
    return(double())
  }
  p <- as.double(p)
  alpha <- as.double(alpha)
  beta <- as.double(beta)
  lambda <- as.double(lambda)
  if (lower.tail) lt_f <- 1L else lt_f <- 0L
  if (log.p) lp_f <- 1L else lp_f <- 0L
  if (exact) {
    QDLAP <- .Call(qdelap_C, p, alpha, beta, lambda, lt_f, lp_f,
                   getDelapThreads())
  } else if (length(alpha) > 1L || length(beta) > 1L || length(lambda) > 1L ||
          anyNA(p)) {
    warning("Quantile approximation relies on pooling and is not accurate when",
            " passed vector-valued parameters, NaNs, or NAs. Using exact",
            " version instead.")
    QDLAP <- .Call(qdelap_C, p, alpha, beta, lambda, lt_f, lp_f,
                   getDelapThreads())
  } else if (alpha <= 0 || beta <= 0 || lambda <= 0) {
    QDLAP <- rep.int(NaN, length(p))
  } else {
      if (log.p) p <- exp(p)
      if (!lower.tail) p <- 1 - p
      validIdx <- p > 0 & p < 1
      QDLAP <- double(length(p))
      QDLAP[p < 0] <- NaN
      QDLAP[p == 0] <- 0
      QDLAP[p >= 1] <- Inf
      if (any(validIdx)) {
        n <- min(10 ^ (ceiling(log(alpha * beta + lambda, 10)) + 5), 1e7)
        shiftedGammas <- rgamma(n, shape = alpha, scale = beta)
        DP <- rpois(n, lambda = (shiftedGammas + lambda))
        QDLAP[validIdx] <- as.vector(quantile(DP, p[validIdx], na.rm = TRUE,
                                              type = 8L))
      }
  }
  if (any(is.nan(QDLAP))) warning("NaNs produced")
  QDLAP
}

rdelap <- function(n, alpha, beta, lambda, exact = TRUE) {
  # Follow base R convention (e.g. rpois): a vector n means length(n)
  # variates, and NA or negative n is an error.
  if (length(n) > 1L) n <- length(n) else n <- as.integer(n)
  if (is.na(n) || n < 0L) {
    stop("invalid arguments")
  }
  # Zero-length parameters produce NAs with a warning, mirroring base R.
  if (length(alpha) == 0L || length(beta) == 0L || length(lambda) == 0L) {
    warning("NaNs produced")
    return(rep.int(NaN, n))
  }
  alpha <- as.double(alpha)
  beta <- as.double(beta)
  lambda <- as.double(lambda)
  if (exact) {
    RDLAP <- .Call(rdelap_C, n, alpha, beta, lambda, getDelapThreads())
  } else if (any(alpha <= 0) || any(beta <= 0) || any(lambda <= 0)) {
    warning("NaNs produced")
    RDLAP <- (rep.int(NaN, n))
  } else {
    shiftedGammas <- rgamma(n, shape = alpha, scale = beta)
    RDLAP <- rpois(n, lambda = (shiftedGammas + lambda))
  }
  if (any(is.nan(RDLAP))) warning("NaNs produced")
  RDLAP
}

MoMdelap <- function(x, type = 2L) { # nolint object_name_linter
  type <- as.integer(type)
  if (!(type %in% c(1L, 2L, 3L))) stop("Skew type must be one of 1, 2, or 3.")
  if (length(x) < 3L) {
    stop("MoMdelap requires at least three data points.")
  }
  moMDLAP <- .Call(MoMdelap_C, as.double(x), type)
  if (any(moMDLAP <= 0)) {
    stop("Method of moments not appropriate for this data; results include ",
         "non-positive parameters.")
  }
  moMDLAP
}
