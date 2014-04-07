\name{NEWS}
\alias{NEWS}

\title{NEWS file for the \code{\link[=Delaporte-package]{Delaporte} package}}

\section{Changes in version 2.1-0 (2014-04-06)}{
  \subsection{New Features}{
    \itemize{
      \item Added function to calculate method of moments estimates for the parameters. These estimates should also serve as good starting points for maximum likelihood estimates.
    }
  }
\section{Changes in version 2.0-1 (2014-02-05)}{
  \subsection{Updates}{
    \itemize{
      \item Minor cleanup after Rcpp update (remove Makevar files; tweak documentation).
    }
  }
}\section{Changes in version 2.0-0 (2014-02-05)}{
  \subsection{Updates}{
    \itemize{
      \item Updating for changes made to Rcpp 0.11.0. Previous versions may no longer compile properly and the current version may no longer work with older versions of Rcpp.
    }
  }
}
\section{Changes in version 1.1-0 (2014-01-17)}{
  \subsection{Updates}{
    \itemize{
      \item Change "Depends" to "LinkingTo" in line with CRAN philosophy on package dependencies
      \item Update documentation to give rough estimate as to when "exact = FALSE" should be used (rarely)
      \item Have R code be written in a consistent style
    }
  }
}
\section{Changes in version 1.0-1 (2013-09-01)}{
  \subsection{Updates}{
    \itemize{
      \item Update README now that CRAN uses it as well
      \item Minor grammar and typography tweaks to documentation
    }
  }
}  
\section{Changes in version 1.0-0 (2013-07-09)}{
  \subsection{New Features}{
    \itemize{
    \item Added Poisson-negative binomial approximation to quantile and random variate functions
      \itemize{
      \item This greatly speeds up calculation for large \eqn{\alpha, \beta}, or \eqn{\lambda}.
      }
    }
  }
}
\section{Changes in version 0.1-2 (2013-04-10)}{
  \subsection{New Features}{
    \itemize{
      \item Converted NEWS to NEWS.Rd
      
      \item Added repository and bug report URL to DESCRIPTION
      
    }
  }  
  \subsection{Bug Fixes}{
    \itemize{  
      \item Help documentation corrections, comments, and tweaks
    }  
  }
}    
\section{Changes in version 0.1-1 (2013-04-08)}{
  \itemize{
    \item Initial release
  }
}