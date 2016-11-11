delp <- function(x, alpha, beta, lambda, lg = FALSE){
  if(!is.double(x)) {storage.mode(x) <- 'double'}
  if(!is.double(alpha)) {storage.mode(alpha) <- 'double'}
  if(!is.double(beta)) {storage.mode(beta) <- 'double'}
  if(!is.double(lambda)) {storage.mode(lambda) <- 'double'}
  .Call('ddelap_f_s_wrap', x, alpha, beta, lambda, lg, PACKAGE = 'DelF')
}

gam_ln <- function(x){
  if(!is.double(x)) {storage.mode(x) <- 'double'}
  .Call('gamln_wrap', x, PACKAGE = 'DelF')
}
