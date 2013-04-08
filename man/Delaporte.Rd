\name{Delaporte}
\alias{Delaporte}
\alias{ddelap}
\alias{pdelap}
\alias{qdelap}
\alias{rdelap}
\title{The Delaporte Distribution}
\description{
Density, distribution function, quantile function and random generation for the negative binomial distribution with parameters \code{alpha}, \code{beta}, and \code{lambda}.
}
\usage{
ddelap(x, alpha, beta, lambda, log = FALSE)
pdelap(q, alpha, beta, lambda, lower.tail = TRUE, log.p = FALSE)
qdelap(p, alpha, beta, lambda, lower.tail = TRUE, log.p = FALSE)
rdelap(n, alpha, beta, lambda)
}
\arguments{
  \item{x}{vector of (non-negative integer) quantiles.}
  \item{q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{alpha}{the alpha parameter of the gamma portion of the Delaporte distribution.
    Must be strictly positive, but need not be integer.}
  \item{beta}{the beta parameter of the gamma portion of the Delaporte distribution.
    Must be strictly positive, but need not be integer.}
  \item{lambda}{the lambda parameter of the Poisson portion of the Delaporte distribution.
    Must be strictly positive, but need not be integer.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
}
\details{
  The Delaporte distribution with parameters \code{alpha}, \code{beta}, and \code{lambda} is a discrete probability distribution which has a probability mass function (called by \code{ddelap}) of:
  \deqn{p(n) = \sum_{i=0}^n\frac{\Gamma(\alpha+i)\beta^i\lambda^{n-i}e^{-\lambda}}{\Gamma(\alpha) i! (1+\beta)^{\alpha+i}(n-i)!}}{p(n) = \sum (i=0:n) [\Gamma(\alpha+i) \beta^i \lambda^(n-i) exp(-\lambda)] / [\Gamma(\alpha) i! (1+\beta)^(\alpha+i) (n-i)!]}
  for \eqn{n = 0, 1, 2, \ldots}; \eqn{\alpha, \beta, \lambda > 0}.
  
  It can be considered the convolution of a negative binomial distribution with a Poisson distribution.  Alternatively, it can be considered a counting distribution with both a Poisson and negative binomial component.
  
  The cumulative distribution function, \code{pdelap}, is calcuated through double summation:
  \deqn{CDF(n) = \sum_{j=0}^n \sum_{i=0}^j\frac{\Gamma(\alpha+i)\beta^i\lambda^{j-i}e^{-\lambda}}{\Gamma(\alpha)i!(1+\beta)^{\alpha+i}(j-i)!}}{CDF(n) = \sum(j=0:n) \sum(i=0:j) [\Gamma(\alpha+i) \beta^i \lambda^(j-i) exp(-\lambda)] / [\Gamma(\alpha) i! (1+\beta)^(\alpha+i) (j-i)!]}
  for \eqn{n = 0, 1, 2, \ldots}; \eqn{\alpha, \beta, \lambda > 0}.  This can be computationally expensive for large \code{alpha}, \code{beta}, and \code{n}.  To optimize for speed, the function identifies the largest value passed to it and computes a vector of CDF values for all integers up to and including that value.  The remaining results are looked up off of this vector, requiring only one double summation instead of \code{length(q)} such summations. 
  
  The quantile function, \code{qdelap}, works by identifying the largest valid quantile passed to it, and then, starting from \code{0}, building a CDF vector until such point as the CDF value exceeds the largest passed quantile.  The function returns \code{NaN} for quantiles \code{< 0}, \code{0} for quantiles \code{= 0}, \code{Inf} for quantiles \code{>= 1}.  For valid quantilies, the function uses the constructed CDF vector, finds the first value where the CDF is greater than the quantile, and returns that value.
  
  The random generator, \code{rdelap}, uses inversion and works very similarly to the quantile function, in that a vector of \code{n} uniformly distributed random variates between \code{0} and \code{1} is generated, a CDF vector is constructed as per the quantile function, and the resulting values are looked up off of the constructed vector.
}
\value{
  \code{ddelap} gives the probability mass function,   \code{pdelap} gives the cumulative distribution function,   \code{qdelap} gives the quantile function, and \code{rdelap} generates random deviates.

  Invalid quantiles passed to \code{qdelap} will result in return values of
  \code{NaN} or \code{Inf} as appropriate.
  
  The length of the result is determined by \code{n} for
  \code{rnbinom}, and is the length of the vector of inputs for the other functions.  
}
\references{
  Johnson, Norman Lloyd, Kemp, Adrienne W. and Kotz, Samuel (2005) \emph{Univariate discrete distributions} (Third ed.). John Wiley & Sons. pp. 241--242. ISBN 978-0-471-27246-5.
  
  Wikipedia contributors, "Delaporte distribution," \emph{Wikipedia, The Free Encyclopedia}, \url{http://en.wikipedia.org/w/index.php?title=Delaporte_distribution&oldid=540316091} (accessed April 8, 2013). 
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \link{Distributions} for standard distributions, including \code{\link{dnbinom}} for the negative binomial distribution
  and \code{\link{dpois}} for the Poisson distribution
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{distribution}

