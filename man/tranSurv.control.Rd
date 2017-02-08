\name{tranSurv.control}
\alias{tranSurv.control}
\title{Auxiliary for Controlling tranSurvfit Fitting}
\description{
  Auxiliary function as user interface for \code{tranSurvfit}
  fitting. 
}
\usage{
tranSurv.control(interval, lower = min(interval), upper = max(interval))
}
\arguments{
  \item{interval}{a vector containing the end-points of the interval to be
    searched the transformation parameter.}
  \item{lower}{the lower end point of the interval to be searched.}
  \item{upper}{the upper end point of the interval to be searched.}
}
%\details{}
\value{
  A list with the arguments as components.
}
\author{
  Sy Han Chiou
}
\seealso{\code{tranSurvfit}.}

