% Copyright (c) 2013, Avraham Adler All rights reserved
% SPDX-License-Identifier: BSD-2-Clause

\name{NEWS}
\alias{NEWS}
\title{NEWS file for the Delaporte package}

\section{Version 8.4.0.9000 (2024-06-30)}{
  \subsection{Changed}{
    \itemize{
      \item Added package metadata tests.
      \item Updated \file{DESCRIPTION}.
      \item Testing for \code{NaN} within Fortran now more efficient.
    }
  }
}

\section{Version 8.4.0 (2024-04-04)}{
  \subsection{Changed}{
    \itemize{
      \item \code{qdelap} now automatically calls the \code{exact = TRUE}
      version when passed inappropriate parameters instead of just stopping with
      an error.
      \item Refactor \file{README.md} into component \acronym{MD} files for
      Github purposes.
      \item Capture and restore original thread count after running tests.
      \item Efficiency tweaks to underlying Fortran code.
    }
  }
}

\section{Version 8.3.0 (2023-10-02)}{
  \subsection{Changed}{
    \itemize{
      \item Converted \code{log1p} implementation to one based on its Taylor
      sequence, instead of catastrophic substitution elimination, as some
      compilers may optimize out that subtraction. Since any difference in
      result should be at machine precision or lower, this is considered a
      minor update.
      \item Removed \code{back} keyword from \code{minloc} call in the Fortran
      implementation of \code{qdelap}. Technically uneeded since its value is
      the default \code{FALSE}. This allows for systems with an only partial
      implementation of Fortran 2008 to compile (see
      \href{https://github.com/aadler/Delaporte/issues/3}{Issue 3} and thanks to
      \href{https://github.com/dpascall}{David Pascall}).
      \item Correct error message of \code{setDelapThreads} to reflect that
      setting threadcount to 0 is also an error. Functionality has not changed.
      \item Removed unnecessary explicit check for \code{NaN} in \code{qdelap}
      as it will be covered by check using \code{anyNA}.
      \item Tweaked documentation for clarity.
    }
  }
  \subsection{Added}{
    \itemize{
      \item Added unit test for new \code{log1p} function.
      \item Added \code{LTO_OPT} and \code{LTO_FC_OPT} flags to \file{Makevars}.
    }
  }
}

\section{Version 8.2.0 (2023-08-08)}{
  \subsection{Added}{
    \itemize{
      \item Added \code{getDelapThreads} and \code{setDelapThreads} to
      explicitly control OpenMP thread usage based on
      \href{https://cran.r-project.org/doc/manuals/R-exts.html#OpenMP-support}{
      WRE 1.2.1.1}.
      \item Added \pkg{parallel} as an import for \code{detectCores}.
      \item Added \acronym{SPDX} header to more files.
      \item Added more unit tests.
    }
  }
  \subsection{Changed}{
    \itemize{
      \item Changed Fortran suffix to \code{.f90} for better compatibility with
      the Intel Fortran compiler, per \acronym{CRAN}.
      \item Tweaked minor linting and documentation issues.
      \item Explicitly enumerated exports instead of using wide-ranging regex.
    }
  }
}

\section{Version 8.1.1 (2023-06-19)}{
  \subsection{Changed}{
    \itemize{
      \item Removed \code{-std=f2008} on the advice of Professor Ripley.
    }
  }
}

\section{Version 8.1.0 (2023-02-01)}{
  \subsection{Changed}{
    \itemize{
      \item Fortran source code now reliant on Fortran 2008.
      \item Fortran module \code{ieee_arithmetic} used.
      \item Updated help documentation.
      \item Updated Github actions.
      \item Updated internal code to be more compliant with personal style guide.
      \item Updated CITATION to use \code{bibentry} per \acronym{CRAN} request.
    }
  }
  \subsection{Removed}{
    \itemize{
      \item \code{gamln} replaced with Fortran 2008 \code{log_gamma}.
      \item \code{position} replaced with Fortran \code{minloc}.
      \item The \code{set_nan} and \code{set_inf} procedures were replaced with
      \code{ieee_arithmetic} calls.
      \item Checks for \code{NaN} now use \code{ieee_arithmetic} calls.
    }
  }
}

\section{Version 8.0.3 (2022-09-05)}{
  \subsection{Changed}{
    \itemize{
      \item Updated package with information necessary to acheive
      \acronym{OpenSSF} silver status.
    }
  }
}

\section{Version 8.0.2 (2021-01-20)}{
  \subsection{Changed}{
    \itemize{
      \item Converted unit testing framework to \pkg{tinytest} to reduce
      dependencies.
      \item Augmented citation with \acronym{DOI} and \acronym{CFF} file.
    }
  }
  \subsection{Fixed}{
    \itemize{
      \item Corrected help file.
    }
  }
}

\section{Version 8.0.1 (2021-01-10)}{
  \subsection{Removed}{
    \itemize{
      \item Removed \acronym{SIMD} calls added in 8.0.0 due to Solaris.
      \item Removed \acronym{LTO} mentions from Makevars added in 8.0.0 per
      Brian Ripley.
    }
  }
  \subsection{Changed}{
    \itemize{
      \item Pass more \R-consistent error messages.
      \item Made more robust to improper inputs such as negative, \code{NA}, or
      \code{NaN} parameters or passed values.
      \item R versions of \code{qdelap} and \code{rdelap} behave more like the
      Fortran versions.
    }
  }
}

\section{Version 8.0.0 (2021-01-05) [YANKED]}{
  \subsection{Changed}{
    \itemize{
      \item Exact code more in line with R defaults for d/p/q/r functions in
      that \code{NaN} inputs return \code{NaN} outputs.
      \item Handling of \code{Inf} and \code{-Inf} now more in line with R
      defaults for d/p/q/r functions.
      \item Parameter errors now check for \eqn{< 0} instead of \code{< EPS}.
      \item \code{ddelap} and \code{pdelap} now use 8-byte integers
      (\acronym{INT64}) for their counters which allows calculation for values
      up to \eqn{2^{63}}{2^63}.
      \itemize{
        \item \code{ddelap} will return 0 for all values \eqn{> 2^{63}}{> 2^63}
        under the assumption that the \acronym{PMF} is spread too thinly.
        \item \code{pdelap} will not run for values \eqn{> 2^{63}}{> 2^63}.
        \item \code{ddelap} of \code{Inf} returns 0 and \code{pdelap} of
        \code{Inf} returns 1.
      }
      \item \code{pdelap} gives a warning for finite values
      \eqn{> 2^{15}}{> 2^15} as, depending on the parameters, this may take
      hours to finish (it's a double summation!). The user has to respond with
      \kbd{y} to continue.
      \item Updated help and tests.
      \item Allow compiler to use OpenMP \acronym{SIMD} constructs.
      \item Allow using \acronym{LTO}.
      \item Switched continuous integration from Travis/Appveyor to Github
      actions (thanks to Jim Hester for his help in fixing the coverage check
      issue).
    }
  }
}

\section{Version 7.0.5 (2020-11-10)}{
  \subsection{Fixed}{
    \itemize{
      \item Corrected \acronym{CRAN} requests.
      \item Cleaned up NEWS.
    }
  }
}

\section{Version 7.0.4 (2020-11-08)}{
  \subsection{Changed}{
    \itemize{
      \item Converted tests to be compliant with \pkg{testthat}: 3rd edition.
      \item Cleaned up prose.
    }
  }
  \subsection{Added}{
    \itemize{
      \item Added a citation to DESCRIPTION.
      \item Added a test to increase code coverage.
    }
  }
}

\section{Version 7.0.3 (2020-06-01)}{
  \subsection{Added}{
    \itemize{
      \item CITATION file.
      \item Explicit citation information to README.
      \item Travis, AppVeyor, and CodeCov integration.
      \item Added \code{R_forceSymbols(dll, TRUE)} to source per \acronym{WRE}
      5.4.
      \item Tests for package versioning.
    }
  }
  \subsection{Changed}{
    \itemize{
      \item Moved canonical repository to Github.
      \item Use \code{FFLAGS} instead of \code{FCFLAGS} in Makevars.
    }
  }
  \subsection{Removed}{
    \itemize{
      \item Removed depsy badge from README as site is now defunct.
    }
  }
}

\section{Version 7.0.2 (2019-05-17)}{
  \subsection{Changed}{
    \itemize{
      \item Updated NEWS to reflect recent changes.
      \item Updated various file versions to reflect recent changes.
    }
  }
}

\section{Version 7.0.1 (2019-04-26)}{
  \subsection{Changed}{
    \itemize{
      \item Makevars adjusted by \acronym{CRAN} to prevent issues with parallel
      make.
    }
  }
}

\section{Version 7.0.0 (2019-02-12)}{
  \subsection{Changed}{
    \itemize{
      \item Changed Makevars flags and targets and forced dependance on \R 3.6+
      to comply with new \acronym{CRAN} policies on Fortran and OpenMP.
      \item Preventing spurious floating point errors by implementing hard floor
      of 0 and hard ceiling of 1.
      \item Updated comments in \file{delaporte.f95} source code.
      \item Changed \code{rdelap} testing to always use a fixed seed.
      \item Documentation cleaned up.
    }
  }
  \subsection{Removed}{
    \itemize{
      \item Removed the \code{cleanzeros} function.
      \item Removed unneeded Makevars target \code{cleanup}.
      \item Removed \acronym{GNU}-specific Makevars target \code{.NOTPARALLEL}.
    }
  }
}

\section{Version 6.3.0 (2018-11-21)}{
  \subsection{Added}{
    \itemize{
      \item Added \code{cleanzeros} function that zaps to 0 any value less than
      or equal to \acronym{EPS} (thanks Kipper Fletez-Brant).
      \item Test for zeros (thanks Kipper Fletez-Brant).
      \item Added items to README to comply with
      \href{https://bestpractices.coreinfrastructure.org/en}{\acronym{CII}}
      best practices and added badge.
      \item Added various badges to README.
      \item Added copyright and license information inside \R source code file.
    }
  }
}

\section{Version 6.2.0 (2018-06-22)}{
  \subsection{Added}{
    \itemize{
      \item \code{MoMdelap} now allows a skew \code{type} selection, similar to
      that of the \link[e1071]{skewness} function in the \pkg{e1071} package
      (thanks to William Cipolli for the suggestion). Skew defaults to former
      selection so there should be no change to any values.
      \item Tests added for all three skewness variants.
      \item \acronym{ORCID} added to author description.
      \item Copyright date corrected to reflect actual copyright and not date of
      licensure change.
      }
  }
}

\section{Version 6.1.0 (2017-08-13)}{
  \subsection{Fixed}{
    \itemize{
      \item There was an error in the calculation of \code{beta}, and thus all
      three parameters, in \code{MoMdelap} (thanks to William Cipolli).
      }
  }
}

\section{Version 6.0.0 (2017-03-31)}{
  \subsection{Fixed}{
    \itemize{
      \item The approximate version of \code{qdelap} was not correct. It pooled
      together ALL generated random variables. As it is impractical to generate
      millions of random variates for each unique parameter triplet, the
      functionality was changed so that only the \dQuote{exact} version of
      \code{qdelap} may be used if any of the parameters are a vector.
    }
  }
  \subsection{Changed}{
    \itemize{
      \item In the continued struggle to have the package pass all its tests on
      Solaris \acronym{SPARC}, all logicals needing to be passed to Fortran are
      now passed as integers and not booleans. While all other platforms had no
      problem using Fortran 2003's ISO_C_bindings, Solaris \acronym{SPARC} did.
    }
  }
  \subsection{Added}{
    \itemize{
      \item Added note to README about Solaris \acronym{SPARC}.
      \item More tests including for singleton approximate \code{qdelap}.
    }
  }  
  \subsection{Removed}{
    \itemize{
      \item Older implementation of \code{qdelap} and \code{rdelap}
      approximations are now completely removed from code and documentation.
    }
  }
}

\section{Version 5.0.1 (2017-03-05)}{
  \subsection{Changed}{
    \itemize{
      \item Solaris, being big-endian, appears unable to use the same bit-mask
      transfer for \code{NaN} and \code{Inf} as any other \R platform.
      Therefore, the package has to revert to calling C to call the \R magic
      words which should recognize the platform.
    }
  }
}

\section{Version 5.0.0 (2017-03-05)}{
  \subsection{Changed}{
    \itemize{
      \item Changed handling of non-integers to be in line with \R convention
      that they have 0 probability instead of being considered the next integer.
      \code{ddelap} will return 0 for each non-integral value and will pass one
      warning, unlike base \R which passes one warning per non-integral value.
      \code{qdelap} will just return the value for the last integer and no
      warnings, similar to base \R.
      \item Use bit-representation of \code{NaN} and \code{Inf} to remove need
      to call C for setting these values. This allows the log-gamma and
      singleton d/p/qdelap functions to be elemental.
      \item Changed registration of the C \acronym{SEXP} functions seen by the
      package to use \code{R_CallMethodDef} as now preferred by \acronym{CRAN}.
      \item Explicitly listed the five C \acronym{SEXP} functions seen by the
      package as \code{R_RegisterCCallable} for potential downstream package
      use.
      \item Corrected some prose in documentation page.
      \item Older NEWS entries refactored for accuracy.
    }
  }
  \subsection{Fixed}{
    \itemize{
      \item Explicitly cast \code{n} as an integer in loggamma routine. This
      should suppress warning on \acronym{OSX} mavericks.
    }
  }
  \subsection{Removed}{
    \itemize{
      \item Older implementation of \code{qdelap} and \code{rdelap}
      approximations is now defunct. Only the new faster approximation based on
      section 3.9 of Karlis & Xekalaki (2005) remains. The next minor version
      may remove the \code{.Defunct} calls and the references to the old version
      in the documentation.
      \item Removed \code{colClasses = 'double'} from tests. This should solve
      Sun \acronym{SPARC} issues (It didn't).
    }
  }
  \subsection{Added}{
    \itemize{
      \item Specified, clarified, and added many more unit tests.
    }
  }
}

\section{Version 4.0.3 (2017-01-30)}{
  \subsection{Changed}{
    \itemize{
      \item Back end code ported from C++ & \pkg{Rcpp} to Fortran and C. In
      most cases this alone showed an increase in speed. There should be no
      noticeable difference in results. Random numbers generated with the same
      seeds should be the same as in previous packages.
      \item OpenMP calls activated for systems which have it which allows
      parallel processing of certain loops resulting in further speed increases
      over and above those supplied by Fortran. OpenMP compatibility is
      approximately version 3.1; no \acronym{SIMD} instructions are used.
      \item License changed from \acronym{GPL2+/LGPL3+} to \acronym{BSD-2}
      clause.
    }
  }
  \subsection{Fixed}{
    \itemize{
      \item Includes various emergency fixes to address parallel make and other
      compilation issues in versions 4.0.0--4.0.2.
    }
  }
}

\section{Version 4.0.2 (2017-01-30) [YANKED]}{
  \itemize{
    \item Did not build properly on some systems using a parallel make.
  }
}

\section{Version 4.0.1 (2017-01-30) [YANKED]}{
  \itemize{
    \item Did not build properly on Fedora using Clang.
  }
}

\section{Version 4.0.0 (2017-01-29) [YANKED]}{
  \itemize{
    \item Did not build properly on various platforms such as Fedora and Solaris.
  }
}

\section{Version 3.0.0 (2016-07-19)}{
  \subsection{Changed}{
    \itemize{
      \item The distributional functions now return \code{NaN} where the passed
      \code{alpha}, \code{beta}, or \code{lambda} is close to 0 (less than
      machine epsilon). Within vectors, valid triplets should still return
      proper values. The approximate versions of \code{qdelap} and \code{rdelap}
      do not handle 0 properly, and an error is thrown directing the user to the
      exact version (thanks are due to Jonathan Cairns and Jin-Rui Xu for their
      help in testing).
      }
    }
  \subsection{Added}{
    \itemize{
      \item Faster versions of the approximate \code{qdelap} and \code{rdelap}
      have been implemented.
      \item Tests for near 0 handling on exact version.
      }
    }
  \subsection{Deprecated}{
    \itemize{
      \item The old approximation versions of \code{qdelap} and \code{rdelap}
      can still be used by passing \code{old = TRUE}. This will be removed
      eventually.
    }
  }
}

\section{Version 2.3.1 (2016-06-16)}{
  \subsection{Changed}{
    \itemize{
      \item Reverted back to 2.2-3 as change to trapping zeros caused serious
      downstream issues. Will work on more elegant correction shortly
      (thanks to Mikhail Spivakov, Jonathan Cairns, and Jin-Rui Xu).
    }
  }
}

\section{Version 2.3.0 (2016-06-02)}{
  \subsection{Deprecated}{
    \itemize{
      \item In the next release (\strong{3.0.0}), passing \code{FALSE} to
      \code{exact} in \code{rdelap} and \code{qdelap} will use the method
      described in section 3.9 of Karlis & Xekalaki (2005) to generate
      approximately Delaporte random variates much more quickly for large
      \code{n} and large values of the parameters. However, when implemented,
      this will break reproducible even for specifically set seeds. The new code
      is commented out in the function.
    }
  }
  \subsection{Changed}{
    \itemize{
      \item Change non-exact \code{qdelap} mode to \emph{numeric} instead of
      \emph{integer} to handle returned \code{Inf}.
      \item Refactor NEWS to more closely comply with
      \href{https://keepachangelog.com/}{keep a changelog} suggestions.
    }
  }
  \subsection{Fixed}{
    \itemize{
      \item Check for any parameters not strictly greater than 0 (thanks to
      Mikhail Spivakov and Jonathan Cairns).
    }
  }
}

\section{Version 2.2-3 (2015-06-30)}{
  \subsection{Fixed}{
    \itemize{
      \item Explicitly import code from package \code{\link{stats}}.
      \item Updated code in unit tests to reflect change in
      \pkg{testthat} package (move to \code{expect_equal}).
      \item Minor tweaks to prose and grammar and expanded description.
    }
  }
}

\section{Version 2.2-2 (2014-12-04)}{
  \subsection{Fixed}{
    \itemize{
      \item Updated DESCRIPTION file to comply with recent changes.
    }
  }
}

\section{Version 2.2-1 (2014-09-24)}{
  \subsection{Fixed}{
    \itemize{
      \item Trapped machine precision issues when calling right-tailed version
      of cumulative distribution function (thanks to Mikhail Spivakov and
      Jonathan Cairns).
      \item Update documentation to reflect rounding procedure when non-integer
      passed as value.
    }
  }
}

\section{Version 2.2-0 (2014-07-15)}{
  \subsection{Added}{
    \itemize{
      \item Vectorized \code{d/p/q/rdelap}; \code{alpha}, \code{beta}, and
      \code{lambda} can now be vectors.
      \item Added more robust tests.
      \item Minor optimizations for speed.
    }
  }
  \subsection{Fixed}{
    \itemize{
      \item Corrected \code{lower.tail} and \code{log.p} which were not
      implemented properly.
    }
  }
}

\section{Version 2.1-1 (2014-04-07)}{
  \subsection{Fixed}{
    \itemize{
      \item Adjusted unit test so it can be run more easily on 32 bit machine
      with limited memory.
      \item Fixed spelling and other prose.
    }
  }
}

\section{Version 2.1-0 (2014-04-06)}{
  \subsection{Added}{
    \itemize{
      \item Added function to calculate method of moments estimates for the
      parameters. These estimates should also serve as good starting points for
      maximum likelihood estimation.
      \item Added unit testing.
    }
  }
}

\section{Version 2.0-1 (2014-02-05)}{
  \subsection{Changed}{
    \itemize{
      \item Minor cleanup after \pkg{Rcpp} update (remove \file{Makevar} files
      and tweak documentation).
    }
  }
}

\section{Version 2.0-0 (2014-02-05)}{
  \subsection{Changed}{
    \itemize{
      \item Updating for changes made to \pkg{Rcpp} 0.11.0. Previous versions
      may no longer compile properly and the current version may no longer work
      with older versions of \pkg{Rcpp}.
    }
  }
}

\section{Version 1.1-0 (2014-01-17)}{
  \subsection{Changed}{
    \itemize{
      \item Change \dQuote{Depends} to \dQuote{LinkingTo} in line with
      \acronym{CRAN} philosophy on package dependencies.
      \item Update documentation to give rough estimate as to when
      \code{exact = FALSE} should be used (rarely).
      \item Have \R code be written in a consistent style.
    }
  }
}

\section{Version 1.0-1 (2013-09-01)}{
  \subsection{Changed}{
    \itemize{
      \item Update README now that \acronym{CRAN} uses it as well.
      \item Minor grammar and typography tweaks to documentation.
    }
  }
}

\section{Version 1.0-0 (2013-07-09)}{
  \subsection{Added}{
    \itemize{
      \item Poisson-negative binomial approximation to quantile and random
      variate functions. This greatly speeds up calculation for large
      \eqn{\alpha, \beta}, or \eqn{\lambda}.
    }
  }
}

\section{Version 0.1-2 (2013-04-10)}{
  \subsection{Added}{
    \itemize{
      \item Converted NEWS to NEWS.Rd.
      \item Added repository and bug report \acronym{URL} to DESCRIPTION.
    }
  }  
  \subsection{Fixed}{
    \itemize{  
      \item Help documentation corrections, comments, and tweaks.
    }  
  }
}

\section{Version 0.1-1 (2013-04-08)}{
  \itemize{
    \item Initial version.
  }
}
