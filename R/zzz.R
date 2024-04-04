# nocov start
.onLoad <- function(libname, pkgname) {
  # Get reasonable selction for "max" cpus
  options("DLPCPU" = as.integer( # nolint
    getOption("DLPCPU", default = parallel::detectCores())
  ))
}

.onDetach <- function(libpath) {
  # Restore reasonable option and remove holding variable
  setDelapThreads(getOption("DLPCPU"))
  options("DLPCPU" = NULL) # nolint
}

.onUnload <- function(libpath) {
  library.dynam.unload("Delaporte", libpath)
}
# nocov end
