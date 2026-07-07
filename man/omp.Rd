% Copyright (c) 2023, Avraham Adler All rights reserved
% SPDX-License-Identifier: BSD-2-Clause

\name{setDelapThreads}
\alias{setDelapThreads}
\alias{getDelapThreads}
\alias{openMP}
\alias{openmp}
\title{
Set or get the number of OpenMP threads \pkg{Delaporte} should use
}
\description{
Set or get number OpenMP threads to be used by \pkg{Delaporte} functions which
are parallelized. These include \code{ddelap}, \code{pdelap}, \code{qdelap}, and
\code{rdelap}.
}
\usage{
setDelapThreads(n)
getDelapThreads()
}
\arguments{
  \item{n}{integer: maximum number of threads to be used.}
}

\details{
The requested thread count is stored in an environment local to the package
namespace, following the guidance in Writing R Extensions, and is passed to the
compiled routines on each call, where it scopes only \pkg{Delaporte}'s own
parallel regions. Neither the global options registry nor the process-wide
OpenMP settings are touched, so other OpenMP-using packages are unaffected.

On load, the thread count is initialized to the OpenMP runtime maximum at that
moment, which respects environment variables such as \env{OMP_NUM_THREADS} and
\env{OMP_THREAD_LIMIT}.

For systems which return valid values for
\code{\link[parallel]{detectCores}}, the maximum number of threads will be
capped at that value.
}
\value{
\code{getDelapThreads} returns an integer representing the maximum number of
allowed threads.

\code{setDelapThreads} invisibly returns the integer value actually set after
flooring and capping.
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
