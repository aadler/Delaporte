ddelap <-
  function (x, alpha, beta, lambda, log = FALSE) 
  {
    DDLAP <- ddelap_C(x, alpha, beta, lambda, log)
    return(DDLAP)
  }
pdelap <-
  function (q, alpha, beta, lambda, lower.tail = TRUE, log.p = FALSE) 
  {
    PDLAP <- pdelap_C(q, alpha, beta, lambda, lower.tail, log.p)
    return(PDLAP)
  }
qdelap <-
  function (p, alpha, beta, lambda, lower.tail = TRUE, log.p = FALSE) 
  {
    QDLAP <- qdelap_C(p, alpha, beta, lambda, lower.tail, log.p)
    return(QDLAP)
  }
rdelap <-
  function (n, alpha, beta, lambda) 
  {
    RDLAP <- rdelap_C(n, alpha, beta, lambda)
    return(RDLAP)
  }