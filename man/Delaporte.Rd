\name{Delaporte}
\alias{Delaporte}
\alias{ddelap}
\alias{pdelap}
\alias{qdelap}
\alias{rdelap}
\title{The Delaporte Distribution}
\description{
Density, distribution function, quantile function and random generation for the negative binomial distribution with parameters \code{alpha}, \code{beta} and \code{lambda}.
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
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
