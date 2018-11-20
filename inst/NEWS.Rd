\name{NEWS}
\alias{NEWS}
\title{NEWS file for the \code{\link[=Delaporte-package]{Delaporte} package}}

\section{Version 7.0.0 (2018-11-20)}{
  \subsection{Added}{
    \itemize{
      \item Added "cleanzeros" function that zaps to 0 any value less than a margin times EPS. That margin is currently 1.1. This does mean that any values formerly between 2.2204460492503131e-16 and 2.4424906541753446e-16 will now return 0. This is a breaking change, thus the major model increment, but it should not affect anyone drastically. See \href{https://bitbucket.org/aadler/delaporte/issues/1/negative-p-values-in-pdelap}{Issue #1} (Thanks Kipper Fletez-Brant).
      \item Test for zeros. Using call submitted in \href{https://bitbucket.org/aadler/delaporte/issues/1/negative-p-values-in-pdelap}{Issue #1} (Thanks Kipper Fletez-Brant).
      \item Added items to READMe to comply with \href{https://bestpractices.coreinfrastructure.org/en}{CII} best practices and added badge.
      \item Added \code{SHLIB_OPENMP_CFLAGS} to \code{PKG_LIBS} to try and address \code{NOTE} on r-devel+linux even though \href{https://cran.r-project.org/doc/manuals/r-release/R-exts.html#OpenMP-support}{Writing R extensions 1.2.1.1} specifically states "This does not apply to Fortran 9x code, where SHLIB_OPENMP_FCFLAGS should be used in both PKG_FCFLAGS and PKG_LIBS."
    }
  }
}

\section{Version 6.2.0 (2018-06-22)}{
  \subsection{Added}{
    \itemize{
      \item \code{MoMdelap} now allows a skew \code{type} selection, similar to that of the \link[e1071]{skewness} function in the \pkg{e1071} package (thanks to William Cipolli for the suggestion). Skew defaults to former selection so there should be no change to any values.
      \item Tests added for all three skewness variants.
      \item ORCID added to author description
      \item Copyright date corrected to reflect actual copyright and not date of licensure change.
      }
  }
}

\section{Version 6.1.0 (2017-08-13)}{
  \subsection{Fixed}{
    \itemize{
      \item There was an error in the calculation of \code{beta}, and thus all three parameters, in \code{MoMdelap} (thanks to William Cipolli).
      }
  }
}

\section{Version 6.0.0 (2017-03-31)}{
  \subsection{Fixed}{
    \itemize{
      \item The approximate version of \code{qdelap} was not correct. It pooled together ALL generated random variables. As it is impractical to generate millions of random variates for each unique parameter triplet, the functionality was changed so that only the \dQuote{exact} version of \code{qdelap} may be used if any of the parameters are a vector.
    }
  }
  \subsection{Changed}{
    \itemize{
      \item In the continued struggle to have the package pass all its tests on Solaris SPARC, all logicals needing to be passed to Fortran are now passed as integers and not booleans. While all other platforms had no problem using Fortran 2003's ISO_C_bindings, Solaris SPARC did.
    }
  }
  \subsection{Added}{
    \itemize{
      \item Added note to README about Solaris SPARC.
      \item More tests including for singleton approximate \code{qdelap}.
    }
  }  
  \subsection{Removed}{
    \itemize{
      \item Older implementation of \code{qdelap} and \code{rdelap} approximations are now completely removed from code and documentation.
    }
  }
}

\section{Version 5.0.1 (2017-03-05)}{
  \subsection{Changed}{
    \itemize{
      \item Solaris, being big-endian, appears unable to use the same bit-mask transfer for NaN and INF as any other R platform. Therefore, the package has to revert to calling C to call the R magic words which should recognize the platform.
    }
  }
}

\section{Version 5.0.0 (2017-03-05)}{
  \subsection{Changed}{
    \itemize{
      \item Changed handling of non-integers to be in line with R convention that they have 0 probability instead of being considered the next integer. \code{ddelap} will return 0 for each non-integral value and will pass one warning, unlike base R which passes one warning per non-integral value. \code{qdelap} will just return the value for the last integer and no warnings, similar to base R.
      \item Use bit-representation of NaN and INF to remove need to call C for setting these values. This allows the log-gamma, and singleton d/p/qdelap functions to be elemental.
      \item Changed registration of the C SEXP functions seen by the package to use \code{R_CallMethodDef} as now preferred by CRAN.
      \item Explicitly listed the five C SEXP functions seen by the package as R_RegisterCCallable for potential downstream package use.
      \item Corrected some prose in documentation page.
      \item Older NEWS entries refactored for accuracy.
    }
  }
  \subsection{Fixed}{
    \itemize{
      \item Explicitly cast \code{n} as an integer in loggamma routine. This should suppress warning on OSX mavericks.
    }
  }
  \subsection{Removed}{
    \itemize{
      \item Older implementation of \code{qdelap} and \code{rdelap} approximations is now defunct. Only the new faster approximation based on section 3.9 of Karlis & Xekalaki (2005) remains. The next minor version may remove the \code{.Defunct} calls and the references to the old version in the documentation.
      \item Removed \code{colClasses = 'double'} from tests. This should solve Sun Sparc issues. (It didn't)
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
      \item Back end code ported from C++ & Rcpp to Fortran and C. In most cases this alone showed an increase in speed. There should be no noticeable difference in results. Random numbers generated with the same seeds should be the same as in previous packages.
      \item OpenMP calls activated for systems which have it which allows parallel processing of certain loops resulting in further speed increases over and above those supplied by Fortran. OpenMP compatibility is approximately version 3.1; no SIMD instructions are used.
      \item License changed from GPL2+/LGPL3+ to BSD-2 clause.
    }
  }
  \subsection{Fixed}{
    \itemize{
      \item Includes various emergency fixes to address parallel make and other compilation issues in versions 4.0.0--4.0.2.
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
    \item Did not build properly on various platforms such as Fedora and Solaris
  }
}

\section{Version 3.0.0 (2016-07-19)}{
  \subsection{Changed}{
    \itemize{
      \item The distributional functions now return \code{NaN} where the passed \code{alpha}, \code{beta}, or \code{lambda} is close to 0 (less than machine epsilon). Within vectors, valid triplets should still return proper values. The approximate versions of \code{qdelap} and \code{rdelap} do not handle 0 properly, and an error is thrown directing the user to the exact version. Thanks are due to Jonathan Cairns and Jin-Rui Xu for their help in testing.
      }
    }
  \subsection{Added}{
    \itemize{
      \item Faster versions of the approximate \code{qdelap} and \code{rdelap} have been implemented.
      \item Tests for near 0 handling on exact version.
      }
    }
  \subsection{Deprecated}{
    \itemize{
      \item The old approximation versions of \code{qdelap} and \code{rdelap} can still be used by passing \code{old = TRUE}. This will be removed eventually.
    }
  }
}

\section{Version 2.3.1 (2016-06-16)}{
  \subsection{Changed}{
    \itemize{
      \item Reverted back to 2.2-3 as change to trapping zeros caused serious downstream issues. Will work on more elegant correction shortly (thanks to Mikhail Spivakov, Jonathan Cairns, and Jin-Rui Xu).
    }
  }
}

\section{Version 2.3.0 (2016-06-02)}{
  \subsection{Deprecated}{
    \itemize{
      \item In the next release (\strong{3.0.0}), passing FALSE to \code{exact} in \code{rdelap} and \code{qdelap} will use the method described in section 3.9 of Karlis & Xekalaki (2005) to generate approximately Delaporte random variates much more quickly for large \code{n} and large values of the parameters. However, when implemented, this will break reproducible even for specifically set seeds. The new code is commented out in the function.
    }
  }
  \subsection{Changed}{
    \itemize{
      \item Change non-exact \code{qdelap} mode to \emph{numeric} instead of \emph{integer} to handle returned \code{Inf}.
      \item Refactor NEWS to more closely comply with \href{http://keepachangelog.com/}{keep a changelog} suggestions.
    }
  }
  \subsection{Fixed}{
    \itemize{
      \item Check for any parameters not strictly greater than 0 (thanks to Mikhail Spivakov and Jonathan Cairns).
    }
  }
}

\section{Version 2.2-3 (2015-06-30)}{
  \subsection{Fixed}{
    \itemize{
      \item Explicitly import code from package \code{\link{stats}}.
      \item Updated code in unit tests to reflect change in \code{\link{testthat}} package (move to \code{expect_equal}).
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
      \item Trapped machine precision issues when calling right-tailed version of cumulative distribution function (thanks to Mikhail Spivakov and Jonathan Cairns).
      \item Update documentation to reflect rounding procedure when non-integer passed as value.
    }
  }
}

\section{Version 2.2-0 (2014-07-15)}{
  \subsection{Added}{
    \itemize{
      \item Vectorized \code{d/p/q/rdelap}; \code{alpha}, \code{beta}, and \code{lambda} can now be vectors.
      \item Added more robust tests.
      \item Minor optimizations for speed.
    }
  }
  \subsection{Fixed}{
    \itemize{
      \item Corrected \code{lower.tail} and \code{log.p} which were not implemented properly.
    }
  }
}

\section{Version 2.1-1 (2014-04-07)}{
  \subsection{Fixed}{
    \itemize{
      \item Adjusted unit test so it can be run more easily on 32 bit machine with limited memory.
      \item Fixed spelling and other prose.
    }
  }
}

\section{Version 2.1-0 (2014-04-06)}{
  \subsection{Added}{
    \itemize{
      \item Added function to calculate method of moments estimates for the parameters. These estimates should also serve as good starting points for maximum likelihood estimation.
      \item Added unit testing.
    }
  }
}

\section{Version 2.0-1 (2014-02-05)}{
  \subsection{Changed}{
    \itemize{
      \item Minor cleanup after Rcpp update (remove Makevar files; tweak documentation).
    }
  }
}

\section{Version 2.0-0 (2014-02-05)}{
  \subsection{Changed}{
    \itemize{
      \item Updating for changes made to Rcpp 0.11.0. Previous versions may no longer compile properly and the current version may no longer work with older versions of Rcpp.
    }
  }
}

\section{Version 1.1-0 (2014-01-17)}{
  \subsection{Changed}{
    \itemize{
      \item Change "Depends" to "LinkingTo" in line with CRAN philosophy on package dependencies.
      \item Update documentation to give rough estimate as to when "exact = FALSE" should be used (rarely).
      \item Have R code be written in a consistent style.
    }
  }
}

\section{Version 1.0-1 (2013-09-01)}{
  \subsection{Changed}{
    \itemize{
      \item Update README now that CRAN uses it as well.
      \item Minor grammar and typography tweaks to documentation.
    }
  }
}

\section{Version 1.0-0 (2013-07-09)}{
  \subsection{Added}{
    \itemize{
      \item Poisson-negative binomial approximation to quantile and random variate functions.
        \itemize{
          \item This greatly speeds up calculation for large \eqn{\alpha, \beta}, or \eqn{\lambda}.
        }
    }
  }
}

\section{Version 0.1-2 (2013-04-10)}{
  \subsection{Added}{
    \itemize{
      \item Converted NEWS to NEWS.Rd.
      \item Added repository and bug report URL to DESCRIPTION.
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