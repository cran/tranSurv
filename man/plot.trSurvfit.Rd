% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.trSurvfit.R
\name{plot.trSurvfit}
\alias{plot.trSurvfit}
\title{Plot the survival estimation based on the structural transformation model}
\usage{
\method{plot}{trSurvfit}(x, ...)
}
\arguments{
\item{x}{an object of class \code{trSurvfit} returned by the \code{trSurvfit()} or the \code{trReg()} function.}

\item{...}{graphical parameters to be passed to methods.}
}
\value{
A \code{ggplot} object.
}
\description{
Plot the survival estimation for an \code{trSurvfit}/\code{trReg} object.
}
\examples{
data(channing, package = "boot")
chan <- subset(channing, entry < exit)

plot(trReg(Surv(entry, exit, cens) ~ 1, data = chan))
}
\keyword{Plot}
