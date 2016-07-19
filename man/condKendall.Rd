\name{condKendall}
\alias{condKendall}
\title{
  Test for quasi-independence on truncated survival data
}
\description{
  This function performs statistical test for quasi-independence
  between truncation and failure time. The hypothesis test is based on
  the conditional Kendall's tau of Martin and Betensky (2005) or two
  versions of the inverse probability weighted Kendall's tau of Austin
  and Betensky (2014). 
}
\usage{
condKendall(trun, obs, delta = NULL, method = "MB", weights = NULL,
a = 0, trans = "linear", ...)
}
\arguments{
  \item{trun}{ left truncation time, satisfying trun <= obs}
  \item{obs}{observed failure time, must be the same length as
    \code{trun}, might be right-censored.}
  \item{delta}{an optional vector of censoring indicator (0 = censored,1
    = event) for obs. When this vector is not specified, the
    function assumes there is no censoring and all observed failure time
    experienced an event.}
  \item{method}{a character string specifying the different version of
    conditional Kendall's tau. The following are permitted:
    \code{MB}: conditional Kendall's tau proposed in Martin and Betensky
    (2005) as \deqn{\hat{\tau_c},}
    \code{IPW1}: inverse probability weighted estimator proposed in
    Austin and Betensky (2014) as \deqn{\hat{\tau_{c2}},}
    \code{IPW2}: restricted inverse probability weighted estimator proposed in
    Austin and Betensky (2014) as \deqn{\hat{\tau_{c3}}.}}
  \item{weights}{ an optional vector of sampling weights used when
    \code{weights=IPW1} or \code{weights=IPW2}. Default weights will be
    computed when leaving it as null.}
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
%% \details{}
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
  Martin E. and Betensky R. A. (2005), Testing quasi-independence of failure
  and truncation times via conditional Kendall's tau, \emph{Journal of the
    American Statistical Association}, \bold{100} (470): 484-492.

  Austin, M. D. and Betensky R. A. (2014), Eliminating bias due to
  censoring in Kendall's tau estimators for quasi-independence of
  truncation and failure, \emph{Computational Statistics & Data
    Analysis}, \bold{73}: 16-26.
}
% \seealso{\code{\link{permTest}}}
\examples{
## Generate simulated data from transformation model
datgen <- function(n) {
    a <- -0.3
    X <- rweibull(n, 2, 4) ## failure times
    U <- rweibull(n, 2, 1) ## latent truncation time
    T <- (1 + a) * U - a * X ## apply transformation
    C <- rlnorm(n, .8, 1) ## censoring
    dat <- data.frame(trun = T, obs = pmin(X, C), delta = 1 * (X <= C))
    return(subset(dat, trun <= obs))
}

set.seed(123)
dat <- datgen(300)

with(dat, condKendall(trun, obs, delta))
with(dat, condKendall(trun, obs, delta, method = "IPW1"))
with(dat, condKendall(trun, obs, delta, method = "IPW2"))
}

