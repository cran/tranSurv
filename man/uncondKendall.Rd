\name{uncondKendall}
\alias{uncondKendall}
\title{
  Computes the unconditional Kendall's tau
}
\description{
  This function computes the unconditional Kendall's tau, or the Kendall
  rank correlation coefficient for two variables. The implementation is
  done in C and performs a lot faster than \code{cor} with \code{method
    = "kendall"}.
}
\usage{
uncondKendall(x, y)
}
\arguments{
  \item{x}{ a numeric vector}
  \item{y}{ a numeric vector with compatible length to 'y'}
}
%% \details{}
\value{
  The output is a numeric value for the unconditional Kendall's tau. 
}
\references{
  Kendall, M. G. (1938), A new measure of rank correlation,
  \emph{Biometrika}, 81-93. 
}
% \seealso{\code{\link{permTest}}}
\examples{
## Prepare to generate from normal copula
pho <- .25
A <- matrix(c(1, pho, pho, sqrt(1 - pho^2)), 2)

## Generating truncated survival data from normal copula
set.seed(123)
n <- 10000
dat <- matrix(NA, n, 2)
for (i in 1:n) dat[i,] <- pnorm(A \%*\% rnorm(2))
dat <- data.frame(dat)
colnames(dat) <- c("x", "y")
system.time(print(uncondKendall(dat$x, dat$y)))
system.time(print(cor(dat$x, dat$y, method = "kendall")))
}

