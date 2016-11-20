# Delaporte

**Delaporte** is an `R` package which provides the probability mass, distribution, quantile, random variate generation, and method of moments parameter estimation functions for the Delaporte distribution. As the distribution does not have a closed form but requires summations or double summations to calculate values, the functions have been programmed in Fortran and C. In cases where approximations are sufficient, the quantile and random variate generator have the option to use a much faster Poisson-negative binomial estimate as opposed to the full Delaporte double summations.

The author is grateful to [Drew Schmidt](https://github.com/wrathematics) both generally for his writings on R, C++, and Fortran and specifically for help with this project.

####Note
 - This project attempts to follow [Semantic Versioning](http://semver.org/)
 - This project attempts to follow the changelog system at [Keep a CHANGELOG](http://keepachangelog.com/)