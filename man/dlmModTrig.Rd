\name{dlmModTrig}
\alias{dlmModTrig}
\title{Create Fourier representation of a periodic DLM component}
\description{
  The function creates a dlm representing a specified periodic
  component. 
}
\usage{
dlmModTrig(s, q, om, tau, dV = 1, dW = 0, m0, C0)
}
\arguments{
  \item{s}{the period, if integer.}
  \item{q}{number of harmonics in the DLM.}
  \item{om}{the frequency.}
  \item{tau}{the period, if not an integer.}
  \item{dV}{variance of the observation noise.}
  \item{dW}{a single number expressing the variance of the system noise.}
  \item{m0}{\eqn{m_0}{m0}, the expected value of the pre-sample state vector.}
  \item{C0}{\eqn{C_0}{C0}, the variance matrix of the pre-sample state vector.}
}
\details{
  The periodic component is specified by one and only one of \code{s},
  \code{om}, and \code{tau}. When \code{s} is given, the function
  assumes that the period is an integer, while a period specified by
  \code{tau} is assumed to be noninteger. Instead of \code{tau},
  the frequency \code{om} can be specified. The argument \code{q}
  specifies the number of harmonics to include in the model. When
  \code{tau} or \code{omega} is given, then \code{q} is required as
  well, since in this case the implied Fourier representation has
  infinitely many harmonics. On the other hand, if \code{s} is given,
  \code{q} defaults to all the harmonics in the Fourier representation,
  that is \code{floor(s/2)}.

  The system variance of the resulting dlm is \code{dW} times the identity
  matrix of the appropriate dimension. 
}
\value{
  An object of class dlm, representing a periodic component.
}
\references{Petris, Petrone, and Campagnoli, Dynamic Linear Models with
  R, Springer (2009). West and Harrison, Bayesian forecasting and dynamic models
  (2nd ed.), Springer, 1997.}
\author{Giovanni Petris \email{GPetris@uark.edu}}

\seealso{\code{\link{dlmModSeas}}, \code{\link{dlmModARMA}},
  \code{\link{dlmModPoly}}, \code{\link{dlmModReg}}} 
\examples{
dlmModTrig(s = 3)
dlmModTrig(tau = 3, q = 1) # same thing
dlmModTrig(s = 4) # for quarterly data
dlmModTrig(s = 4, q = 1)
dlmModTrig(tau = 4, q = 2) # a bad idea!
m1 <- dlmModTrig(tau = 6.3, q = 2); m1
m2 <- dlmModTrig(om = 2 * pi / 6.3, q = 2)
all.equal(unlist(m1), unlist(m2))
}
\keyword{misc}% at least one, from doc/KEYWORDS

