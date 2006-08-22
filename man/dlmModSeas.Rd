\name{dlmModSeas}
\alias{dlmModSeas}
\title{Create a DLM for seasonal factors}
\description{
  The function creates a DLM representation of seasonal component. 
}
\usage{
dlmModSeas(frequency, m0 = rep(0, frequency - 1), C0 = 1e+07 * diag(nrow
= frequency - 1), dV = 1, dW = c(1, rep(0, frequency - 2))) 
}
\arguments{
  \item{frequency}{How many seasons?}
  \item{m0}{\eqn{m0}.}
  \item{C0}{\eqn{C0}.}
  \item{dV}{Variance of the observational noise.}
  \item{dW}{Diagonal elements of the variance matrix of the system noise.}
}
\value{
  An object of class dlm representing a seasonal factor for a process
  with \code{frequency} seasons.
}
\references{Harvey, Forecasting, structural time series models and the
  Kalman filter, Cambridge University Press, 1989.}
\author{Giovanni Petris, \email{GPetris@uark.edu}}

\seealso{\code{\link{dlmModARMA}}, \code{\link{dlmModPoly}},
  \code{\link{dlmModReg}}, and \code{\link{dlmModTrig}} for the Fourier
  representation of a seasonal component.}
\examples{
## seasonal component for quarterly data
dlmModSeas(4, dV = 3.2)
}
\keyword{misc}% at least one, from doc/KEYWORDS

