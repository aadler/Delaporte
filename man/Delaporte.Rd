\name{Delaporte}
\alias{Delaporte}
\alias{delaporte}
\alias{ddelap}
\alias{pdelap}
\alias{qdelap}
\alias{rdelap}
\alias{MoMdelap}
\title{The Delaporte Distribution}
\description{
Density, distribution, quantile, random variate generation, and method of moments parameter estimation functions for the Delaporte distribution with parameters \code{alpha}, \code{beta}, and \code{lambda}. 
}
\usage{
ddelap(x, alpha, beta, lambda, log = FALSE)
pdelap(q, alpha, beta, lambda, lower.tail = TRUE, log.p = FALSE)
qdelap(p, alpha, beta, lambda, lower.tail = TRUE, log.p = FALSE, exact=TRUE)
rdelap(n, alpha, beta, lambda, exact = TRUE)

MoMdelap(x, type = 2L)
}
\arguments{
  \item{x}{vector of (non-negative integer) quantiles.}
  \item{q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.}
  \item{alpha}{vector of alpha parameters of the gamma portion of the Delaporte
  distribution.
    Must be strictly positive, but need not be integer. }
  \item{beta}{vector of beta parameters of the gamma portion of the Delaporte
  distribution.
    Must be strictly positive, but need not be integer.}
  \item{lambda}{vector of lambda parameters of the Poisson portion of the
  Delaporte distribution.
    Must be strictly positive, but need not be integer.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
  \item{exact}{logical; if TRUE uses double summation to generate quantiles or
  random variates. Otherwise uses Poisson-negative binomial approximation.}
  \item{type}{integer; 1L will return \strong{g1}, 2L will return \strong{G1},
  and 3L will return \strong{b1}, as per \code{\link[e1071]{skewness}}.}
}
\details{
The Delaporte distribution with parameters \code{alpha}, \code{beta}, and
\code{lambda} is a discrete probability distribution which can be considered the
convolution of a negative binomial distribution with a Poisson distribution.
Alternatively, it can be considered a counting distribution with both Poisson
and negative binomial components. The Delaporte's probability mass function,
called via \code{ddelap}, is:

  \deqn{p(n) = \sum_{i=0}^n\frac{\Gamma(\alpha+i)\beta^i\lambda^{n-i}e^{-\lambda}}{\Gamma(\alpha) i! (1+\beta)^{\alpha+i}(n-i)!}}{p(n) = \sum (i=0:n) [\Gamma(\alpha+i) \beta^i \lambda^(n-i) exp(-\lambda)] / [\Gamma(\alpha) i! (1+\beta)^(\alpha+i) (n-i)!]}

for \eqn{n = 0, 1, 2, \ldots}; \eqn{\alpha, \beta, \lambda > 0}.
  
 If an element of \code{x} is not integer, the result of \code{ddelap} is zero
with a warning.
  
The Delaporte's cumulative distribution function, \code{pdelap}, is calculated
through double summation:

  \deqn{CDF(n) = \sum_{j=0}^n \sum_{i=0}^j\frac{\Gamma(\alpha+i)\beta^i\lambda^{j-i}e^{-\lambda}}{\Gamma(\alpha)i!(1+\beta)^{\alpha+i}(j-i)!}}{CDF(n) = \sum(j=0:n) \sum(i=0:j) [\Gamma(\alpha+i) \beta^i \lambda^(j-i) exp(-\lambda)] / [\Gamma(\alpha) i! (1+\beta)^(\alpha+i) (j-i)!]}

for \eqn{n = 0, 1, 2, \ldots}; \eqn{\alpha, \beta, \lambda > 0}. If only
singleton values for the parameters are passed in, the function uses a shortcut.
It identifies the largest value passed to it, computes a vector of CDF values
for all integers up to and including that value, and reads the remaining results
from this vector. This requires only one double summation instead of
\code{length(q)} such summations. If at least one of the parameters is itself a
vector of length greater than 1, the function has to build the double summation
for each entry in \code{q}.
  
The quantile function, \code{qdelap}, is right continuous:
\code{qdelap(q, alpha, beta, lambda)} is the smallest integer \eqn{x} such that
\eqn{P(X \le x) \ge q}. This has function has two versions. When
\code{exact = TRUE}, the function builds a CDF vector and the first value for
which the CDF is greater than or equal to \code{q} is returned as the quantile.
While this procedure is accurate, for sufficiently large \eqn{\alpha, \beta},
or \eqn{\lambda} it can take a very long time. Therefore, when dealing with
singleton parameters, \code{exact = FALSE} can be passed to take advantage of
the Delaporte's definition as a counting distribution with both a Poisson and a
negative binomial component. Based on Karlis & Xekalaki (2005) it will generate
\code{n} gamma variates \eqn{\Gamma} with shape \eqn{\alpha} and scale
\eqn{\beta} and then \code{n} pseudo-Delaporte variates as Poisson random
variables with parameter \eqn{\lambda + \Gamma}, finally calling the
\code{\link{quantile}} function on the result. The \dQuote{exact} method is
always more accurate and is also significantly faster for reasonable values of
the parameters. Also, the \dQuote{exact} method must be used when passing
parameter vectors, as the pooling would become intractable. Ad-hoc testing
indicates that the \dQuote{exact} method should always be used until
\eqn{\alpha\beta + \lambda \approx 2500}{\alpha\beta + \lambda ~ 2500}. Both
versions return \code{NaN} for quantiles \eqn{< 0}, \eqn{0} for quantiles
\eqn{= 0}, and \code{Inf} for quantiles \eqn{\ge 1}.
  
The random variate generator, \code{rdelap}, also has multiple versions. When
\code{exact = TRUE}, it uses inversion by creating a vector of \code{n}
uniformly distributed random variates between \eqn{0} and \eqn{1}. If all the
parameters are singletons, a single CDF vector is constructed as per the
quantile function, and the entries corresponding to the uniform variates are
read off of the constructed vector. If the parameters are themselves vectors,
then it passes the entire uniform variate vector to \code{qdelap}, which is
slower. When \code{exact = FALSE}, regardless of the length of the parameters,
it generates \code{n} gamma variates \eqn{\Gamma} with shape \eqn{\alpha} and
scale \eqn{\beta} and then \code{n} pseudo-Delaporte variates as Poisson random
variables with parameter \eqn{\lambda + \Gamma}. As there is no pooling, each
individual random variate reflects the parameter triplet which generated it. The
non-inversion method is usually faster.

\code{MoMdelap} uses the definition of the Delaporte's mean, variance, and skew
to calculate the method of moments estimates of \eqn{\alpha}, \eqn{\beta}, and
\eqn{\lambda}, which it returns as a numeric vector. This estimate is also a
reasonable starting point for maximum likelihood estimation using nonlinear
optimizers such as \code{\link{optim}} or \code{\link[nloptr]{nloptr}}. If the
data is clustered near 0, there are times when method of moments would result in
a non-positive parameter. In these cases \code{MoMdelap} will throw an error.
For the sample skew, the user has the choice to select \strong{g1}, \strong{G1},
or \strong{b1} as defined in Joanes & Gill (1997) and found in
\code{\link[e1071]{skewness}}. The selection defaults to option 2, \strong{G1},
which Joanes & Gill found to have the least mean-square error for non-normal
distributions.
}
\value{
  \code{ddelap} gives the probability mass function, \code{pdelap} gives the
  cumulative distribution function, \code{qdelap} gives the quantile function,
  and \code{rdelap} generates random deviates. Values close to 0 (less than or
  equal to machine epsilon) for \eqn{\alpha, \beta} or \eqn{\lambda} will return
  \code{NaN} for that particular entry. Proper triplets within a set of vectors
  should still return valid values. For the approximate versions of
  \code{qdelap} and \code{rdelap}, having a value close to 0 will trip an error,
  sending the user to the exact version which currently properly handles
  vector-based inputs which contain 0.

  Invalid quantiles passed to \code{qdelap} will result in return values of
  \code{NaN} or \code{Inf} as appropriate.
  
  The length of the result is determined by \code{x} for \code{ddelap}, \code{q}
  for \code{pdelap}, \code{p} for \code{qdelap}, and \code{n} for \code{rdelap}.
  The distributional parameters (\eqn{\alpha, \beta, \lambda}) are recycled as
  necessary to the length of the result.
  
  When using the \code{lower.tail = FALSE} or \code{log / log.p = TRUE} options,
  some accuracy may be lost at knot points or the tail ends of the distributions
  due to the limitations of floating point representation.
  
  \code{MoMdelap} returns a triplet comprising a method-of-moments based
  estimate of \eqn{\alpha}, \eqn{\beta}, and \eqn{\lambda}.
}
\author{
Avraham Adler \email{Avraham.Adler@gmail.com}
}
\references{

  Joanes, D. N. and Gill, C. A. (1998) Comparing Measures of Sample Skewness and
  Kurtosis. \emph{Journal of the Royal Statistical Society. Series D
  (The Statistician)} \bold{47}(1), 183--189.
  \url{https://www.jstor.org/stable/2988433}
  
  Johnson, N. L., Kemp, A. W. and Kotz, S. (2005)
  \emph{Univariate discrete distributions} (Third ed.). John Wiley & Sons.
  pp. 241--242. ISBN 978-0-471-27246-5.
  
  Karlis, D. and Xekalaki, E. (2005) Mixed Poisson Distributions.
  \emph{International Statistical Review} \bold{73}(1), 35--58.
  \url{https://projecteuclid.org/euclid.isr/1112304811}

  Vose, D. (2008) \emph{Risk analysis: a quantitative guide} (Third, illustrated
  ed.). John Wiley & Sons. pp. 618--619 ISBN 978-0-470-51284-5
}  

\seealso{
  \link{Distributions} for standard distributions, including
  \code{\link{dnbinom}} for the negative binomial distribution and
  \code{\link{dpois}} for the Poisson distribution, and
  \code{\link[e1071]{skewness}} for skew options.
}

\examples{
## Density and distribution
A <- c(0, seq_len(50))
PMF <- ddelap(A, alpha = 3, beta = 4, lambda = 10)
CDF <- pdelap(A, alpha = 3, beta = 4, lambda = 10)

##Quantile
A <- seq(0,.95, .05)
qdelap(A, alpha = 3, beta = 4, lambda = 10)
A <- c(-1, A, 1, 2)
qdelap(A, alpha = 3, beta = 4, lambda = 10)

## Compare a Poisson, negative binomial, and three Delaporte distributions with the same mean:
P <- rpois(25000, 25)                                  ## Will have the tightest spread
DP1 <- rdelap(10000, alpha = 2, beta = 2, lambda = 21) ## Close to the Poisson
DP2 <- rdelap(10000, alpha = 3, beta = 4, lambda = 13) ## In between
DP3 <- rdelap(10000, alpha = 4, beta = 5, lambda = 5)  ## Close to the Negative Binomial
NB <- rnbinom(10000, size = 5, mu = 25)                ## Will have the widest spread
mean(P);mean(NB);mean(DP1);mean(DP2);mean(DP3)         ## Means should all be near 25
MoMdelap(DP1);MoMdelap(DP2);MoMdelap(DP3)              ## Estimates should be close to originals
\dontshow{
plot(density(P), col = "black", lwd = 2, main = "Distribution Comparison", xlab = "Value", xlim = c(0, 80))
lines(density(DP1), col = "blue", lwd = 2)
lines(density(DP2), col = "green3", lwd = 2)
lines(density(DP3), col = "orange3", lwd = 2)
lines(density(NB), col = "red", lwd = 2)
legend(x = "topright", legend = c("Poisson {l=25}", "DP {a=2, b=2, l=21}", "DP {a=3, b=4, l=13}",
  "DP {a=4, b=5, l=5}", "NegBinom {a=5, b=5}"), col=c("black", "blue", "green3","orange3", "red"), lwd = 2)
  }
}
\keyword{distribution}