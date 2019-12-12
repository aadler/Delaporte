[![](https://www.r-pkg.org/badges/version-last-release/Delaporte)]
[![](https://cranlogs.r-pkg.org/badges/Delaporte)]
[![Research software impact](http://depsy.org/api/package/cran/Delaporte/badge.svg)](http://depsy.org/package/r/Delaporte)
[![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/2011/badge)](https://bestpractices.coreinfrastructure.org/projects/2011)

# Delaporte

**Delaporte** is an `R` package which provides the probability mass,
distribution, quantile, random variate generation, and method of moments
parameter estimation functions for the Delaporte distribution. As the
distribution does not have a closed form but requires summations or double
summations to calculate values, the functions have been programmed in Fortran
and C. In cases where approximations are sufficient, the quantile and random
variate generator have the option to use a much faster Poisson-negative binomial
estimate as opposed to the full Delaporte double summations.

## Citation
If you use the package, please cite it as:

  Avraham Adler (2013). Delaporte: Statistical Functions for the Delaporte
  Distribution. R package version 7.0.3.
  https://CRAN.R-project.org/package=Delaporte

A BibTeX entry for LaTeX users is:

```
  @Manual{,
    title = {Delaporte: Statistical Functions for the Delaporte Distribution},
    author = {Avraham Adler},
    year = {2013},
    note = {R package version 7.0.3},
    url = {https://CRAN.R-project.org/package=Delaporte},
  }
```

## Acknowledgements
The author is grateful to [Drew Schmidt](https://github.com/wrathematics) both
generally for his writings on R, C++, and Fortran and specifically for help with
this project.

## Contributions
Please ensure that all contributions comply with both [R and CRAN standards for packages](https://cran.r-project.org/doc/manuals/r-release/R-exts.html).

### Versioning
This project attempts to follow [Semantic Versioning](http://semver.org/)

### Changelog
This project attempts to follow the changelog system at
[Keep a CHANGELOG](http://keepachangelog.com/)

### Dependancies
This project intends to have as few dependancies as possible. Please consider
that when writing code.

### Style
Please review and conform to the current code stylistic choices (e.g. 80
character lines, two-space indentations).

### Documentation
Please provide valid .Rd files and **not** roxygen-style documentation.

### Tests
Please review the current test suite and supply similar `testthat`-compatible
unit tests for all added functionality.

### Submission
If you would like to contribute to the project, it may be prudent to first
contact the maintainer via email. A request or suggestion may be raised as an
issue as well. To supply a pull request (PR), please:

 1. Fork the project into your own local Bitbucket repository
 2. Create a branch in your repository in which you will make your changes
 3. Push that branch to Bitbucket and then create a pull request
 
At this point, the PR will be discussed and eventually accepted or rejected.