ddelap <-
  function (x, alpha, beta, lambda, log = FALSE) {
    DDLAP <- vector(length = length(x), mode = "numeric")
    DDLAP <- ddelap_C(x, alpha, beta, lambda, log)
    return(DDLAP)
  }
pdelap <-
  function (q, alpha, beta, lambda, lower.tail = TRUE, log.p = FALSE) {
    PDLAP <- vector(length = length(q), mode = "numeric")
    PDLAP <- pdelap_C(q, alpha, beta, lambda, lower.tail, log.p)
    return(PDLAP)
  }
qdelap <-
  function (p, alpha, beta, lambda, lower.tail = TRUE, log.p = FALSE, exact = TRUE) {
    QDLAP <- vector(length = length(p), mode = "numeric")
    if (exact) {
      QDLAP <- qdelap_C(p, alpha, beta, lambda, lower.tail, log.p)
    } else {
      pValid <- p[p > 0 & p < 1]
      pNan <- p[p < 0]
      p0 <- p[p == 0]
      pInf <- p[p >= 1]
      n <- min(10 ^ (ceiling(log(alpha * beta + lambda, 10)) + 3), 1e7)
      NB <- rnbinom(n, mu = alpha * beta, size = alpha)
      P <- rpois(n, lambda = lambda)
      DP <- NB + P
      QValid <- as.vector(quantile(DP, pValid, na.rm = TRUE))
      QNan <- rep.int(NaN, times = length(pNan))
      Q0 <- rep.int(0, times = length(p0))
      QInf <- rep.int(Inf, times = length(pInf))
      QDLAP <- as.vector(c(QNan, Q0, QValid, QInf))
    }  
    return(QDLAP)
  }
rdelap <-
  function (n, alpha, beta, lambda, exact = TRUE) {
    RDLAP <- vector(length = length(n), mode = "numeric")
    if (exact) {
      RDLAP <- rdelap_C(n, alpha, beta, lambda)
    } else {
      NB <- rnbinom(max(1e7, n), mu = alpha * beta, size = alpha)
      P <- rpois(max(1e7, n), lambda = lambda)
      DP <- NB + P
      if (n > 1e7) {
        RDLAP <- DP
      } else {
        RDLAP <- sample(x = DP, size = n, replace = TRUE)
      }
    }
    return(RDLAP)
  }

MoMdelap <- function (x) {
    MoMDLAP <- vector(length = 3, mode = "numeric")
    MoMDLAP <- MoMdelap_C(x)
    if (any(MoMDLAP < 0)) stop ("Method of moments not appropriate for this data; results include negative parameters.")
    return(MoMDLAP)
  }