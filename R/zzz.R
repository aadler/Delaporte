# nocov start
.onLoad <- function(libname, pkgname) {
  # Get reasonable selection for "max" cpus and store in enviornment (per 88406)
  DelaporteEnv <<- new.env(parent = emptyenv())
  if (!exists("DLPCPU", envir = DelaporteEnv)) {
    assign("DLPCPU", parallel::detectCores(), envir = DelaporteEnv)
  }
}

.onDetach <- function(libpath) {
  # Restore reasonable option and remove holding variable and environment.
  setDelapThreads(get("DLPCPU", envir = DelaporteEnv))
}

.onUnload <- function(libpath) {
  rm(DelaporteEnv, pos = 1)
  library.dynam.unload("Delaporte", libpath)
}
# nocov end
