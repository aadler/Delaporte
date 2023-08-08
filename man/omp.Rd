% Copyright (c) 2023, Avraham Adler All rights reserved
% SPDX-License-Identifier: BSD-2-Clause

\name{getDelapThreads}
\alias{getDelapThreads}
\alias{setDelapThreads}
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
getDelapThreads()
setDelapThreads(n)
}
\arguments{
  \item{n}{integer: maximum number of threads to be used.}
}

\details{
For systems which return valid values for \code{\link[parallel]{detectCores}},
the maximum number of threads will be capped at that value.
}
\value{
\code{getDelapThreads} returns an integer representing the maximum number of
allowed threads.

\code{setDelapThreads} invisibly returns the integer passed to it.
}
\author{Avraham Adler \email{Avraham.Adler@gmail.com}}
\seealso{
 \code{\link[parallel]{detectCores}} in \pkg{parallel} package.
}
\examples{
getDelapThreads()
setDelapThreads(2L)
getDelapThreads()
}

\keyword{data}
