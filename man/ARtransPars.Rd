\name{ARtransPars}
\alias{ARtransPars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to parametrize a stationary AR process}
\description{%
  The function maps a vector of length p to the vector of autoregressive
  coefficients of a stationary AR(p) process. It can be used to
  parametrize a stationary AR(p) process 
}
\usage{
ARtransPars(raw)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{raw}{a vector of length p}
}
\details{
  The function first maps each element of \code{raw} to (0,1) using
  tanh. The numbers obtained are treated as the first partial
  autocorrelations of a stationary AR(p) process and the vector of the
  corresponding autoregressive coefficients is computed and returned.
}
\value{
  The vector of autoregressive coefficients of a stationary AR(p) process
  corresponding to the parameters in \code{raw}. 
}
\references{Jones, 1987. Randomly choosing parameters from the
  stationarity and invertibility region of autoregressive-moving average
  models. \emph{Applied Statistics}, 36.}
\author{Giovanni Petris, \email{GPetris@uark.edu}}

\examples{
(ar <- ARtransPars(rnorm(5)))
all( Mod(polyroot(c(1,-ar))) > 1 ) # TRUE
}

\keyword{misc}

