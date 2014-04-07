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
Density, distribution, quantile, random variate generation, and method of momenst parameter estimation functions for the Delaporte distribution with parameters \code{alpha}, \code{beta}, and \code{lambda}. 
}
\usage{
ddelap(x, alpha, beta, lambda, log = FALSE)
pdelap(q, alpha, beta, lambda, lower.tail = TRUE, log.p = FALSE)
qdelap(p, alpha, beta, lambda, lower.tail = TRUE, log.p = FALSE, exact = TRUE)
rdelap(n, alpha, beta, lambda, exact = TRUE)

MoMdelap(x)
}
\arguments{
  \item{x}{vector of (non-negative integer) quantiles.}
  \item{q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.}
  \item{alpha}{the alpha parameter of the gamma portion of the Delaporte distribution.
    Must be strictly positive, but need not be integer.}
  \item{beta}{the beta parameter of the gamma portion of the Delaporte distribution.
    Must be strictly positive, but need not be integer.}
  \item{lambda}{the lambda parameter of the Poisson portion of the Delaporte distribution.
    Must be strictly positive, but need not be integer.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are
    \eqn{P[X \le x]}, otherwise, \eqn{P[X > x]}.}
  \item{exact}{logical; if TRUE uses double summation to generate quantiles or random variates. Otherwise uses Poisson-negative binomial approximation.}
}
\details{
  The Delaporte distribution with parameters \code{alpha}, \code{beta}, and \code{lambda} is a discrete probability distribution which can be considered the convolution of a negative binomial distribution with a Poisson distribution.  Alternatively, it can be considered a counting distribution with both Poisson and negative binomial components. The Delaporte's probability mass function, called via \code{ddelap}, is:
  \deqn{p(n) = \sum_{i=0}^n\frac{\Gamma(\alpha+i)\beta^i\lambda^{n-i}e^{-\lambda}}{\Gamma(\alpha) i! (1+\beta)^{\alpha+i}(n-i)!}}{p(n) = \sum (i=0:n) [\Gamma(\alpha+i) \beta^i \lambda^(n-i) exp(-\lambda)] / [\Gamma(\alpha) i! (1+\beta)^(\alpha+i) (n-i)!]}
  for \eqn{n = 0, 1, 2, \ldots}; \eqn{\alpha, \beta, \lambda > 0}.
  
Its cumulative distribution function, \code{pdelap}, is calculated through double summation:
  \deqn{CDF(n) = \sum_{j=0}^n \sum_{i=0}^j\frac{\Gamma(\alpha+i)\beta^i\lambda^{j-i}e^{-\lambda}}{\Gamma(\alpha)i!(1+\beta)^{\alpha+i}(j-i)!}}{CDF(n) = \sum(j=0:n) \sum(i=0:j) [\Gamma(\alpha+i) \beta^i \lambda^(j-i) exp(-\lambda)] / [\Gamma(\alpha) i! (1+\beta)^(\alpha+i) (j-i)!]}
  for \eqn{n = 0, 1, 2, \ldots}; \eqn{\alpha, \beta, \lambda > 0}. This can be computationally expensive for large \code{alpha}, \code{beta}, and \code{n}. Therefore the function identifies the largest value passed to it and computes a vector of CDF values for all integers up to and including that value. The remaining results are read from this vector, requiring only one double summation instead of \code{length(q)} such summations. 
  
The quantile function, \code{qdelap}, has two versions. When \code{exact=TRUE}, the function identifies the largest valid quantile passed to it, and, starting from \eqn{0}, builds a CDF vector until such point as the CDF value exceeds the largest passed quantile. For valid quantiles, the function uses the constructed CDF vector to find the first value for which the CDF is greater than or equal to the passed quantile and returns that value.

However, the above method can take a very long time for large \eqn{\alpha, \beta}, or \eqn{\lambda}. Therefore, when \code{exact=FALSE}, the function takes advantage of the Delaporte's definition as a counting distribution with both a Poisson and negative binomial component. It generates up to \eqn{10^7} variates from a negative binomial distribution with shape \eqn{\alpha} and scale \eqn{\beta}, and up to another \eqn{10^7} variates from a Poisson distribution with the mean \eqn{\lambda}. It then sums the two sets of variates and calls the \code{\link{quantile}} function on the result. Ad-hoc testing indicates that the "exact" option is still significantly faster until \eqn{\alpha\beta + \lambda \approx 5000}{\alpha\beta + \lambda ~ 5000}. Both versions return \code{NaN} for quantiles \eqn{< 0}, \eqn{0} for quantiles \eqn{= 0}, and \code{Inf} for quantiles \eqn{\ge 1}.
  
The random variate generator, \code{rdelap}, has the same two versions as the quantile function. When \code{exact=TRUE}, it uses inversion and works very similarly to the quantile function, in that a vector of \code{n} uniformly distributed random variates between \eqn{0} and \eqn{1} is generated, a CDF vector is constructed as per the quantile function, and the resulting values are read off of the constructed vector. When \code{exact=FALSE}, the larger of \code{n} or \eqn{10^7} variates from both a Poisson and negative binomial with the appropriate parameters are generated and summed. If \code{n} \eqn{< 10^7}, sampling with replacment is used to generate the \code{n} samples from the pool of \eqn{10^7} pseudo-Delaporte variates.

\code{MoMdelap} uses the definition of the Delaporte's mean, variance, and skew to calculate the method of moments estimates of \eqn{\alpha}, \eqn{\beta}, and \eqn{\lambda}, which it returns as a numeric vector. This estimate is also a reasonable starting point for maximum likelihood estimation using nonlinear optimizers such as \code{\link{optim}} or \code{\link[nloptr]{nloptr}}.
}
\value{
  \code{ddelap} gives the probability mass function, \code{pdelap} gives the cumulative distribution function, \code{qdelap} gives the quantile function, and \code{rdelap} generates random deviates.

  Invalid quantiles passed to \code{qdelap} will result in return values of
  \code{NaN} or \code{Inf} as appropriate.
  
  The length of the result is determined by \code{n} for
  \code{rdelap}, and is the length of the vector of inputs for the other functions.  
}
\author{
Avraham Adler \code{<Avraham.Adler@gmail.com>}
}
\references{
  Johnson, Norman Lloyd, Kemp, Adrienne W. and Kotz, Samuel (2005) \emph{Univariate discrete distributions} (Third ed.). John Wiley & Sons. pp. 241--242. ISBN 978-0-471-27246-5.
  
  Vose, David (2008) \emph{Risk analysis: a quantitative guide} (Third, illustrated ed.). John Wiley & Sons. pp. 618--619 ISBN 978-0-470-51284-5
}

\seealso{
  \link{Distributions} for standard distributions, including \code{\link{dnbinom}} for the negative binomial distribution
  and \code{\link{dpois}} for the Poisson distribution.
}

\examples{
## Density and distribution
A <- c(0, seq_len(50))
PMF <- ddelap(A, alpha = 3, beta = 4, lambda = 10)
CDF <- pdelap(A, alpha = 3, beta = 4, lambda = 10)
summary(1 - cumsum(PMF) / CDF) #Relative Error

##Quantile
A <- seq(0,.95, .05)
qdelap(A, alpha = 3, beta = 4, lambda = 10)
qdelap(A, alpha = 3, beta = 4, lambda = 10, exact = FALSE)
A <- c(-1, A, 1, 2)
qdelap(A, alpha = 3, beta = 4, lambda = 10)
qdelap(A, alpha = 3, beta = 4, lambda = 10, exact = FALSE)

## Compare a Poisson, negative binomial, and three Delaporte distributions with the same mean:
P <- rpois(25000, 25)                                  ## Will have the tightest spread
DP1 <- rdelap(25000, alpha = 2, beta = 2, lambda = 21) ## Close to the Poisson
DP2 <- rdelap(25000, alpha = 3, beta = 4, lambda = 13) ## In between
DP3 <- rdelap(25000, alpha = 4, beta = 5, lambda = 5)  ## Close to the Negative Binomial
NB <- rnbinom(25000, size = 5, mu = 25)                ## Will have the widest spread
mean(P);mean(NB);mean(DP1);mean(DP2);mean(DP3)         ## Means should all be near 25
MoMdelap(DP1);MoMdelap(DP2);MoMdelap(DP3)              ## Estimates should be close to originals
\dontshow{
plot(density(P), col = "black", lwd = 2, main = "Distribution Comparison", xlab = "Value")
lines(density(DP1), col = "blue", lwd = 2)
lines(density(DP2), col = "green3", lwd = 2)
lines(density(DP3), col = "orange3", lwd = 2)
lines(density(NB), col = "red", lwd = 2)
legend(x = "topright", legend = c("Poisson {l=25}", "DP {a=2, b=2, l=21}", "DP {a=3, b=4, l=13}",
  "DP {a=4, b=5, l=5}", "NegBinom {a=5, b=5}"), col=c("black", "blue", "green3","orange3", "red"), lwd = 2)
  }
}
\keyword{distribution}