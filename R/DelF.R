ddelp <- function(x, alpha, beta, lambda, lg = FALSE){
  if(!is.double(x)) {storage.mode(x) <- 'double'}
  if(!is.double(alpha)) {storage.mode(alpha) <- 'double'}
  if(!is.double(beta)) {storage.mode(beta) <- 'double'}
  if(!is.double(lambda)) {storage.mode(lambda) <- 'double'}
  .Call('ddelap_f_wrap', x, alpha, beta, lambda, lg, PACKAGE = 'DelF')
}

pdelp <- function(q, alpha, beta, lambda, lt = TRUE, lg = FALSE){
  if(!is.double(q)) {storage.mode(q) <- 'double'}
  if(!is.double(alpha)) {storage.mode(alpha) <- 'double'}
  if(!is.double(beta)) {storage.mode(beta) <- 'double'}
  if(!is.double(lambda)) {storage.mode(lambda) <- 'double'}
  .Call('pdelap_f_wrap', q, alpha, beta, lambda, lt, lg, PACKAGE = 'DelF')
}

qdelp <- function(p, alpha, beta, lambda){
  if(!is.double(p)) {storage.mode(p) <- 'double'}
  if(!is.double(alpha)) {storage.mode(alpha) <- 'double'}
  if(!is.double(beta)) {storage.mode(beta) <- 'double'}
  if(!is.double(lambda)) {storage.mode(lambda) <- 'double'}
  .Call('qdelap_f_wrap', p, alpha, beta, lambda, PACKAGE = 'DelF')
}
