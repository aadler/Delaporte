# Delaporte

**Delaporte** is an `R` package which provides the probability mass, distribution, quantile, random variate generation, and method of moments parameter estimation functions for the Delaporte distribution. As the distribution does not have a closed form but requires summations or double summations to calculate values, the functions have been programmed in Fortran and C. In cases where approximations are sufficient, the quantile and random variate generator have the option to use a much faster Poisson-negative binomial estimate as opposed to the full Delaporte double summations.

The author is grateful to [Drew Schmidt](https://github.com/wrathematics) both generally for his writings on R, C++, and Fortran and specifically for help with this project.

###Solaris SPARC issues
For those using Solaris SPARC machines, please check the [CRAN check results](https://cran.r-project.org/web/checks/check_results_Delaporte.html). While the package may install without complaints, due to SPARC's being big-endian and the idiosyncracies of its Fortran compiler, spurious results _may_ occur. If the most recent checks show ERROR, the only recourse is to use version 3.0.0---the last known version for which Solaris SPARC did not have any known issues. However, please check the NEWS file to see what needed adjustments may be necessary.

####Note
 - This project attempts to follow [Semantic Versioning](http://semver.org/)
 - This project attempts to follow the changelog system at [Keep a CHANGELOG](http://keepachangelog.com/)