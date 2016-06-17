\name{NEWS}
\alias{NEWS}
\title{NEWS file for the \code{\link[=Delaporte-package]{Delaporte} package}}

\section{Changes in version 2.3.1 (2016-06-16)}{
  \subsection{Changed}{
    \itemize{
      \item Reverted back to 2.2-3 as change to trapping zeros caused serious downstream issues. Will work on more elegant correction shortly (thanks to Mikhail Spivakov, Jonathan Cairns, and Jin-Rui Xu)
    }
  }
}

\section{Changes in version 2.3.0 (2016-06-02)}{
  \subsection{Deprecated}{
    \itemize{
      \item In the next release (\strong{3.0.0}), passing FALSE to \code{exact} in \code{rdelap} and \code{qdelap} will use the method described in section 3.9 of Karlis & Xekalaki (2005) to generate approximately Delaporte random variates much more quickly for large \code{n} and large values of the parameters. However, when implemented, this will break reproduceability even for specifically set seeds. The new code is commented out in the function.
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
      \item Check for any parameters not strictly greater than 0 (thanks to Mikhail Spivakov and Jonathan Cairns)
    }
  }
}

\section{Version 2.2-3 (2015-06-30)}{
  \subsection{Fixed}{
    \itemize{
      \item Explicitly import code from package \code{\link{stats}}
      \item Updated code in unit tests to reflect change in \code{\link{testthat}} package (move to \code{expect_equal})
      \item Minor tweaks to prose and grammar and expanded description
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

\section{Changes in version 2.2-0 (2014-07-15)}{
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
      \item Change "Depends" to "LinkingTo" in line with CRAN philosophy on package dependencies
      \item Update documentation to give rough estimate as to when "exact = FALSE" should be used (rarely)
      \item Have R code be written in a consistent style
    }
  }
}

\section{Version 1.0-1 (2013-09-01)}{
  \subsection{Changed}{
    \itemize{
      \item Update README now that CRAN uses it as well
      \item Minor grammar and typography tweaks to documentation
    }
  }
}

\section{Version 1.0-0 (2013-07-09)}{
  \subsection{Added}{
    \itemize{
      \item Poisson-negative binomial approximation to quantile and random variate functions
        \itemize{
          \item This greatly speeds up calculation for large \eqn{\alpha, \beta}, or \eqn{\lambda}.
        }
    }
  }
}

\section{Version 0.1-2 (2013-04-10)}{
  \subsection{Added}{
    \itemize{
      \item Converted NEWS to NEWS.Rd
      \item Added repository and bug report URL to DESCRIPTION
    }
  }  
  \subsection{Fixed}{
    \itemize{  
      \item Help documentation corrections, comments, and tweaks
    }  
  }
}

\section{Version 0.1-1 (2013-04-08)}{
  \itemize{
    \item Initial version
  }
}