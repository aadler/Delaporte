% Copyright (c) 2023, Avraham Adler All rights reserved
% SPDX-License-Identifier: BSD-2-Clause

\name{getThreads}
\alias{getThreads}
\alias{setThreads}
\alias{openMP}
\alias{openmp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set or get number of OpenMP threads \pkg{Delaporte} should use
}
\description{
Set and get number OpenMP threads to be used in \pkg{Delaporte} functions which
are parallelized. These include \code{ddelap}, \code{pdelap}, \code{qdelap}, and
\code{rdelap}.
}
\usage{
getThreads()
setThreads(n)
}
\arguments{
  \item{n}{integer: maximum number of threads to be used.}
}

\details{
For systems which return valid values for \code{\link[parallel]{detectCores}},
the maximum number of threads will be capped at that value.
}
\value{No return value; used for side effects.}
\author{Avraham Adler \email{Avraham.Adler@gmail.com}}
\seealso{
 \code{\link[parallel]{detectCores}} in \pkg{parallel} package.
}
\examples{
getThreads()
setThreads(2L)
getThreads()
}

\keyword{data}
