# Delaporte #

**Delaporte** is an R package to provide the probability mass, distribution, and quantile functions, as well as a random number generator, for the Delaporte distribution. As the distribution is not closed and requires summations, or, at times, double summations, to calculate values, the functions have been programmed in C++ using the [Rcpp](http://cran.r-project.org/web/packages/Rcpp/index.html) package written by [Dirk Eddelbuettel](http://dirk.eddelbuettel.com/code/rcpp.html) among others.