# Copyright (c) 2013, Avraham Adler All rights reserved
# SPDX-License-Identifier: BSD-2-Clause

plt0err <- 'Parameters must be strictly greater than 0.'

ddelap <- function(x, alpha, beta, lambda, log = FALSE) {
  x <- as.double(x)
  if (any(x[!is.nan(x)] > floor(x[!is.nan(x)]))) {
    warning("Non-integers passed to ddelap. These will have 0 probability.")
  }
  if (log) log_f <- 1L else log_f <- 0L
  .Call(ddelap_C, x, as.double(alpha), as.double(beta), as.double(lambda), log_f)
}

pdelap <- function(q, alpha, beta, lambda, lower.tail = TRUE, log.p = FALSE) {
  if (lower.tail) lt_f <- 1L else lt_f <- 0L
  if (log.p) lp_f <- 1L else lp_f <- 0L
  .Call(pdelap_C, as.double(q), as.double(alpha), as.double(beta),
        as.double(lambda), lt_f, lp_f)
}

qdelap <- function(p, alpha, beta, lambda, lower.tail = TRUE, log.p = FALSE, exact = TRUE) {
  p <- as.double(p)
  alpha <- as.double(alpha)
  beta <- as.double(beta)
  lambda <- as.double(lambda)
  QDLAP <- double(length(p))
  if (exact) {
    if (lower.tail) lt_f <- 1L else lt_f <- 0L
    if (log.p) lp_f <- 1L else lp_f <- 0L
    QDLAP <- .Call(qdelap_C, p, alpha, beta, lambda, lt_f, lp_f)
  } else {
    if (any(alpha <= 0) || any(beta <= 0) || any(lambda <= 0))
      stop(plt0err)
    if (length(alpha) > 1 || length(beta) > 1 || length(lambda) > 1)
      stop(paste('Quantile approximation relies on pooling and thus is not',
                 'accurate when passing vector-valued parameters.',
                 'Please use exact version.'))
    if (log.p) p <- exp(p)
    if (!lower.tail) p <- 1 - p
    pValid <- p[p > 0 & p < 1]
    pNan <- p[p < 0]
    p0 <- p[p == 0]
    pInf <- p[p >= 1]
    n <- min(10 ^ (ceiling(log(alpha * beta + lambda, 10)) + 5), 1e7)
    ShiftedGammas <- rgamma(n, shape = alpha, scale = beta)
    DP <- rpois(n, lambda = (ShiftedGammas + lambda))
    QValid <- as.vector(quantile(DP, pValid, na.rm = TRUE, type = 8))
    QNan <- rep.int(NaN, times = length(pNan))
    Q0 <- rep.int(0, times = length(p0))
    QInf <- rep.int(Inf, times = length(pInf))
    QDLAP <- as.vector(c(QNan, Q0, QValid, QInf), mode = 'integer')
  }
  return(QDLAP)
}

rdelap <- function(n, alpha, beta, lambda, exact = TRUE) {
  n <- as.integer(n)
  alpha <- as.double(alpha)
  beta <- as.double(beta)
  lambda <- as.double(lambda)
  RDLAP <- double(length(n))
  if (!exact) {
    if (any(alpha <= 0) || any(beta <= 0) || any(lambda <= 0)) stop(plt0err)
    ShiftedGammas <- rgamma(n, shape = alpha, scale = beta)
    RDLAP <- rpois(n, lambda = (ShiftedGammas + lambda))
  } else {
    RDLAP <- .Call(rdelap_C, n, alpha, beta, lambda)
  }
  return(RDLAP)
}

MoMdelap <- function(x, type = 2L) {
  type <- as.integer(type)
  if (!(type %in% c(1L, 2L, 3L))) stop('Skew type must be one of 1, 2, or 3.')
  MoMDLAP <- double(3)
  MoMDLAP <- .Call(MoMdelap_C, as.double(x), type)
  if (any(MoMDLAP <= 0)) stop("Method of moments not appropriate for this data; results include non-positive parameters.")
  return(MoMDLAP)
}

.onUnload <- function(libpath) {
  library.dynam.unload("Delaporte", libpath) # nocov
}
