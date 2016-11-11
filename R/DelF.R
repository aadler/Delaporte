delp <- function(x, alpha, beta, lambda, lg = FALSE){
  if(!is.double(x)) {storage.mode(x) <- 'double'}
  if(!is.double(alpha)) {storage.mode(alpha) <- 'double'}
  if(!is.double(beta)) {storage.mode(beta) <- 'double'}
  if(!is.double(lambda)) {storage.mode(lambda) <- 'double'}
  .Call('ddelap_f_wrap', x, alpha, beta, lambda, lg, PACKAGE = 'DelF')
}

gam_ln <- function(x){
  if(!is.double(x)) {storage.mode(x) <- 'double'}
  .Call('gamln_wrap', x, PACKAGE = 'DelF')
}

ev <- function(x, y){
  if(!is.double(x)) {storage.mode(x) <- 'double'}
  if(!is.double(y)) {storage.mode(x) <- 'double'}
  .Call('ev_f_wrap', x, y, PACKAGE = 'DelF')
}
