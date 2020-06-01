# Copyright (c) 2013, Avraham Adler All rights reserved
# SPDX-License-Identifier: BSD-2-Clause

ddelap <- function(x, alpha, beta, lambda, log = FALSE){
  if (!is.double(x)) {storage.mode(x) <- 'double'}
  if (!is.double(alpha)) {storage.mode(alpha) <- 'double'}
  if (!is.double(beta)) {storage.mode(beta) <- 'double'}
  if (!is.double(lambda)) {storage.mode(lambda) <- 'double'}
  if (any(x > floor(x))) {
    warning("Non-integers passed to ddelap. These will have 0 probability.")
  }
  if (log) log_f <- 1L else log_f <- 0L
  .Call(ddelap_C, x, alpha, beta, lambda, log_f)
}

pdelap <- function(q, alpha, beta, lambda, lower.tail = TRUE, log.p = FALSE){
  if (!is.double(q)) {storage.mode(q) <- 'double'}
  if (!is.double(alpha)) {storage.mode(alpha) <- 'double'}
  if (!is.double(beta)) {storage.mode(beta) <- 'double'}
  if (!is.double(lambda)) {storage.mode(lambda) <- 'double'}
  if (lower.tail) lt_f <- 1L else lt_f <- 0L
  if (log.p) lp_f <- 1L else lp_f <- 0L
  .Call(pdelap_C, q, alpha, beta, lambda, lt_f, lp_f)
}

qdelap <- function(p, alpha, beta, lambda, lower.tail = TRUE, log.p = FALSE, 
                   exact = TRUE){
  if (!is.double(p)) {storage.mode(p) <- 'double'}
  if (!is.double(alpha)) {storage.mode(alpha) <- 'double'}
  if (!is.double(beta)) {storage.mode(beta) <- 'double'}
  if (!is.double(lambda)) {storage.mode(lambda) <- 'double'}
  QDLAP <- double(length(p))
  if (exact) {
    if (lower.tail) lt_f <- 1L else lt_f <- 0L
    if (log.p) lp_f <- 1L else lp_f <- 0L
    QDLAP <- .Call(qdelap_C, p, alpha, beta, lambda, lt_f, lp_f)
  } else {
    if (any(alpha <= 0) || any(beta <= 0) || any(lambda <= 0))
      stop('Parameters must be strictly greater than 0. Please use exact version, if necessary, to prevent spurious results.')
    if (length(alpha) > 1 || length(beta) > 1 || length(lambda) > 1)
      stop('Quantile approximation relies on pooling and thus is not accurate when passing vector-valued parameters. Please use exact version.')
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

rdelap <- function(n, alpha, beta, lambda, exact = TRUE){
  if (!is.integer(n)) {storage.mode(n) <- 'integer'}
  if (!is.double(alpha)) {storage.mode(alpha) <- 'double'}
  if (!is.double(beta)) {storage.mode(beta) <- 'double'}
  if (!is.double(lambda)) {storage.mode(lambda) <- 'double'}
  RDLAP <- double(length(n))
  if (!exact) {
    if (any(alpha <= 0) || any(beta <= 0) || any(lambda <= 0))
      stop('Parameters must be strictly greater than 0. Please use exact version, if necessary, to prevent spurious results.')
    ShiftedGammas <- rgamma(n, shape = alpha, scale = beta)
    RDLAP <- rpois(n, lambda = (ShiftedGammas + lambda))
  } else {
    RDLAP <- .Call(rdelap_C, n, alpha, beta, lambda)
  }
  return(RDLAP)
}

MoMdelap <- function(x, type = 2L){
  if (!is.double(x)) {storage.mode(x) <- 'double'}
  if (!is.integer(type)) {storage.mode(type) <- 'integer'}
  if (!(type %in% c(1L, 2L, 3L)))
    stop('Skew type must be one of 1, 2, or 3.')
  MoMDLAP <- double(3)
  MoMDLAP <- .Call(MoMdelap_C, x, type)
  if (any(MoMDLAP <= 0)) stop("Method of moments not appropriate for this data; results include non-positive parameters.")
  return(MoMDLAP)
}

.onUnload <- function(libpath) {
  library.dynam.unload("Delaporte", libpath)
}
