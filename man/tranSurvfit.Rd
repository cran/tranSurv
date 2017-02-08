\name{tranSurvfit}
\alias{tranSurvfit}
\title{This function creates survival curves under dependent truncation
  and independent censoring via a structural transformation model.
}
\description{
  A structural transformation model for a latent, quasi-independent
  truncation time as a function of the observed dependent truncation
  time and the event time, and an unknown dependence parameter. 
  The dependence parameter is chosen to either minimize the absolute
  value of the restricted inverse probability weighted Kendall's tau or
  maximize the corresponding p-value. The marginal distribution for the
  truncation time and the event time are completely left unspecified.
}
\usage{
tranSurvfit(trun, obs, delta = NULL, trans = "linear", plots = FALSE,
control = tranSurv.control(), ...) 
}
\arguments{
  \item{trun}{left truncation time, satisfying trun <= obs.}
  \item{obs}{observed failure time, must be the same length as
    \code{trun}, might be right-censored.}
  \item{delta}{an optional vector of censoring indicator (0 = censored,1
    = event) for obs. When this vector is not specified, the
    function assumes there is no censoring and all observed failure time
    experienced an event.}
  \item{trans}{a character string specifying the transformation
    structure. The following are permitted:
    \code{linear}: linear transformation structure, 
    \code{log}: log-linear transformation structure,
    \code{exp}: exponential transformation structure.}
  \item{plots}{an optional logical value; if TRUE, a series of
    diagnostic plots as well as the survival curve for the observed
    failure time will be plotted.}
  \item{control}{controls lower and upper bounds when
    \code{trans} is an user specified function.}
  \item{...}{ for future methods.}
}
\details{
  The structure of the transformation model is of the form:
  \deqn{h(U) = (1 + a)^-1 * (h(T) + ah(X)),}
  where T is the truncation time, X is the observed failure time, U is
  the transformed truncation time that is quasi-independent from X and
  h(.) is a monotonic transformation function. 
  The condition, T < X, is assumed to be satisfied.
  The quasi-independent truncation time, U, is obtained by inverting the
  test for quasi-independence by either minimizing the absolute
  value of the restricted inverse probability weighted Kendall's tau or
  maximize the corresponding p-value.

  At the current version, three transformation structures can be
  specified. \code{trans = "linear"} corresponds to h(X) = 1;
  \code{trans = "log"} corresponds to h(X) = log(X); \code{trans =
  "exp"} corresponds to h(X) = exp(X). 
}
\value{
  The output contains the following components:
  \item{Sy}{ estiamted survival function at the (ordered) observed points.}
  \item{byTau}{ a list contains the estimator of transformation
    parameter. The following are the components:
    \code{par}: the best set of transformation parameter found.
    \code{obj}: the value of the inverse probability weighted Kendall's
    tau corresponding to 'par'. 
  }
  \item{byP}{ a list contains the estimator of transformation
    parameter. The following are the components:
    \code{par}: the best set of transformation parameter found.
    \code{obj}: the value of p-value based on the inverse probability
    weighted Kendall's tau corresponding to 'par'. 
  }
  \item{qind}{ a data frame consists of two quasi-independent variables:
    \code{trun}: the transformed truncation time.
    \code{obs}: the corresponding uncensored failure time. }
}

\references{
  Martin E. and Betensky R. A. (2005), Testing quasi-independence of failure
  and truncation times via conditional Kendall's tau, \emph{Journal of the
    American Statistical Association}, \bold{100} (470): 484-492.
  
  Austin, M. D. and Betensky R. A. (2014), Eliminating bias due to
  censoring in Kendall's tau estimators for quasi-independence of
  truncation and failure, \emph{Computational Statistics & Data
    Analysis}, \bold{73}: 16-26.

  Chiou, S., Austin, M., Qian, J. and Betensky R. A. (2016),
  Transformation model estimation of survival under dependent truncation
  and independent censoring, an unpublished manuscript. 
}

% \seealso{\code{\link{permTest}}}
\examples{
## Generate simulated data from a transformation model
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
fit <- with(dat, tranSurvfit(trun, obs, delta))
fit

## Checking the transformation parameter
fit$byTau$par
fit$byTau$obj
with(dat, condKendall(trun, obs, delta, method = "IPW2", a = fit$byTau$par))$PE

fit$byP$par
fit$byP$obj
with(dat, condKendall(trun, obs, delta, method = "IPW2", a = fit$byP$par))$p.value
}

