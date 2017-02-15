ddelap <- function(x, alpha, beta, lambda, log = FALSE){
  if(!is.double(x)) {storage.mode(x) <- 'double'}
  if(!is.double(alpha)) {storage.mode(alpha) <- 'double'}
  if(!is.double(beta)) {storage.mode(beta) <- 'double'}
  if(!is.double(lambda)) {storage.mode(lambda) <- 'double'}
  if(any(x > floor(x))) {
    warning("Non-integers passed to ddelap. These will have 0 probability.")
  }
  .Call(ddelap_C, x, alpha, beta, lambda, log)
}

pdelap <- function(q, alpha, beta, lambda, lower.tail = TRUE, log.p = FALSE){
  if(!is.double(q)) {storage.mode(q) <- 'double'}
  if(!is.double(alpha)) {storage.mode(alpha) <- 'double'}
  if(!is.double(beta)) {storage.mode(beta) <- 'double'}
  if(!is.double(lambda)) {storage.mode(lambda) <- 'double'}
  .Call(pdelap_C, q, alpha, beta, lambda, lower.tail, log.p)
}

qdelap <- function(p, alpha, beta, lambda, lower.tail = TRUE, log.p = FALSE, exact = TRUE, old = FALSE){
  if(!is.double(p)) {storage.mode(p) <- 'double'}
  if(!is.double(alpha)) {storage.mode(alpha) <- 'double'}
  if(!is.double(beta)) {storage.mode(beta) <- 'double'}
  if(!is.double(lambda)) {storage.mode(lambda) <- 'double'}
  QDLAP <- double(length(p))
  if (exact) {
    QDLAP <- .Call(qdelap_C, p, alpha, beta, lambda, lower.tail, log.p)
  } else {
    if(any(alpha <= 0) || any(beta <= 0) || any(lambda <= 0))
      stop('Parameters must be strictly greater than 0. Please use exact version, if necessary, to prevent spurious results')
    if (log.p) p <- exp(p)
    if (!lower.tail) p <- 1 - p
    pValid <- p[p > 0 & p < 1]
    pNan <- p[p < 0]
    p0 <- p[p == 0]
    pInf <- p[p >= 1]
    n <- min(10 ^ (ceiling(log(alpha * beta + lambda, 10)) + 5), 1e7)
    if (old) {
      .Defunct(msg = 'This option is defunct. Use old = FALSE. The "old" option may be removed at any time and exact = FALSE will default to the new method.')
    } else {
      ShiftedGammas <- rgamma(n, shape = alpha, scale = beta)
      DP <- rpois(n, lambda = (ShiftedGammas + lambda))
    }
    QValid <- as.vector(quantile(DP, pValid, na.rm = TRUE, type = 8))
    QNan <- rep.int(NaN, times = length(pNan))
    Q0 <- rep.int(0, times = length(p0))
    QInf <- rep.int(Inf, times = length(pInf))
    QDLAP <- as.vector(c(QNan, Q0, QValid, QInf), mode = 'integer')
  }
  return(QDLAP)
}

rdelap <- function(n, alpha, beta, lambda, exact = TRUE, old = FALSE){
  if(!is.integer(n)) {storage.mode(n) <- 'integer'}
  if(!is.double(alpha)) {storage.mode(alpha) <- 'double'}
  if(!is.double(beta)) {storage.mode(beta) <- 'double'}
  if(!is.double(lambda)) {storage.mode(lambda) <- 'double'}
  RDLAP <- double(length(n))
  if (!exact) {
    if(any(alpha <= 0) || any(beta <= 0) || any(lambda <= 0))
      stop('Parameters must be strictly greater than 0. Please use exact version, if necessary, to prevent spurious results')
    if (old) {
      .Defunct(msg = 'This option is defunct. Use old = FALSE. The "old" option may be removed at any time and exact = FALSE will default to the new method.')
    } else {
      ShiftedGammas <- rgamma(n, shape = alpha, scale = beta)
      RDLAP <- rpois(n, lambda = (ShiftedGammas + lambda))
    }
  } else {
    RDLAP <- .Call(rdelap_C, n, alpha, beta, lambda)
  }
  return(RDLAP)
}

MoMdelap <- function(x){
  if(!is.double(x)) {storage.mode(x) <- 'double'}
  MoMDLAP <- double(3)
  MoMDLAP <- .Call(MoMdelap_C, x)
  if (any(MoMDLAP <= 0)) stop ("Method of moments not appropriate for this data; results include non-positive parameters.")
  return(MoMDLAP)
}
