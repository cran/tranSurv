\name{pmcc}
\alias{pmcc}
\title{
  Conditional product-moment correlation coefficient
}
\description{
  This function computes the conditional product-moment correlation
  coefficient proposed by Chen et al. (1996). The conditional
  product-moment correlation coefficient uses only the uncensored
  events. 
}
\usage{
pmcc(trun, obs, a = 0, trans = "linear", ...)
}
\arguments{
  \item{trun}{ left truncation time, satisfying trun <= obs}
  \item{obs}{observed failure time, must be the same length as
    \code{trun}, might be right-censored.}
  \item{a}{ a numeric transformation parameter. The default value is 0,
    which applies no transformation. This parameter must be greater than
    -1. See \code{?tranSurvfit} for the transformation model structure.}
  \item{trans}{a character string specifying the transformation
    structure. The following are permitted:
    \code{linear}: linear transformation structure, 
    \code{log}: log-linear transformation structure,
    \code{exp}: exponential transformation structure.}
  \item{...}{ for future methods.}
}
\value{
  The output contains the following components:
  \item{PE}{ consistent point estimate of conditional Kendall's tau.}
  \item{SE}{ asymptotic standard error the conditional Kendall's tau
    estimator.}
  \item{STAT}{ the value of the normal test statistic.}
  \item{p.value}{ the p-value for the test.}
  \item{trans}{ the transformation model.}
  \item{a}{ transformation parameter.}
}
\references{
  Chen, Chen-Hsin and Tsai, Wei-Yann and Chao, Wei-Hsiung (1996), The
  product-moment correlation coefficient and linear regression for
  truncated data, \emph{Journal of the American Statistical
    Association}, \bold{91}436: 1181-1186.
}
\examples{
## Generate simulated data from transformation model
datgen <- function(n) {
    a <- -0.3
    X <- rweibull(n, 2, 4) ## failure times
    U <- rweibull(n, 2, 1) ## latent truncation time
    T <- (1 + a) * U - a * X ## apply transformation
    C <- Inf ## assuming no censoring
    dat <- data.frame(trun = T, obs = pmin(X, C), delta = 1 * (X <= C))
    return(subset(dat, trun <= obs))
}

set.seed(123)
dat <- datgen(300)

with(dat, pmcc(trun, obs))
}