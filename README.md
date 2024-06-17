---
title: Package Delaporte
---

<!-- badges: start -->
[![](https://www.r-pkg.org/badges/version-last-release/Delaporte)](https://cran.r-project.org/package=Delaporte)
[![](http://cranlogs.r-pkg.org/badges/last-month/Delaporte)](https://cran.r-project.org/package=Delaporte)
[![](https://cranlogs.r-pkg.org/badges/grand-total/Delaporte)](https://cran.r-project.org/package=Delaporte)
[![R-CMD-check](https://github.com/aadler/Delaporte/workflows/R-CMD-check/badge.svg)](https://github.com/aadler/Delaporte/actions)
[![Codecov test coverage](https://codecov.io/gh/aadler/Delaporte/branch/master/graph/badge.svg)](https://app.codecov.io/gh/aadler/Delaporte?branch=master)
[![OpenSSF Best Practices](https://bestpractices.coreinfrastructure.org/projects/2011/badge)](https://bestpractices.coreinfrastructure.org/projects/2011)
<!-- badges: end -->

## Description
**Delaporte** is an `R` package which provides the probability mass,
distribution, quantile, random variate generation, and method of moments
parameter estimation functions for the Delaporte distribution. As the
distribution does not have a closed form, but requires summations or double
summations to calculate values, the functions have been programmed in Fortran
and C. In cases where approximations are sufficient, the quantile and random
variate generator have the option to use a much faster Poisson-negative binomial
estimate as opposed to the full Delaporte double summations.

## Acknowledgment
The author is grateful to [Drew Schmidt](https://github.com/wrathematics) both
generally for his writings on R, C++, and Fortran and specifically for help with
this project.

## Citation
If you use the package, please cite it as per
[CITATION](https://CRAN.R-project.org/package=Delaporte/citation.html).

## Contributions
Please see
[CONTRIBUTING](https://github.com/aadler/delaporte/blob/master/CONTRIBUTING.md).

## Security
Please see [SECURITY](https://github.com/aadler/delaporte/blob/master/SECURITY.md).

## Roadmap
### Major

 * There are no plans for major changes at current.
 
### Minor
 
 * There are no plans for minor changes at current.
 