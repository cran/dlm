\name{dlmModSeas}
\alias{dlmModSeas}
\title{Create a DLM for seasonal factors}
\description{
  The function creates a DLM representation of seasonal component. 
}
\usage{
dlmModSeas(frequency, dV = 1, dW = c(1, rep(0, frequency - 2)),
           m0 = rep(0, frequency - 1),
           C0 = 1e+07 * diag(nrow = frequency - 1)) 
}
\arguments{
  \item{frequency}{how many seasons?}
  \item{dV}{variance of the observation noise.}
  \item{dW}{diagonal elements of the variance matrix of the system noise.}
  \item{m0}{\eqn{m_0}{m0}, the expected value of the pre-sample state vector.}
  \item{C0}{\eqn{C_0}{C0}, the variance matrix of the pre-sample state vector.}
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

