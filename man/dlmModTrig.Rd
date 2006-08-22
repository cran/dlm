\name{dlmModTrig}
\alias{dlmModTrig}
\title{Create Fourier representation of a periodic DLM component}
\description{
  The function creates a dlm representing a specified periodic
  component.} 
}
\usage{
dlmModTrig(s, q, om, tau, m0, C0, dV = 1, dW = 0)
}
\arguments{
  \item{s}{The period, if integer.}
  \item{q}{Number of harmonics in the DLM.}
  \item{om}{The frequency.}
  \item{tau}{The period, if not an integer.}
  \item{m0}{\eqn{m_0}{m0}}
  \item{C0}{\eqn{C_0}{C0}}
  \item{dV}{Variance of the observational noise.}
  \item{dW}{A single number expressing the variance of the system
    noise. }
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
\references{West and Harrison, Bayesian forecasting and dynamic models
  (2nd ed.), Springer, 1997.}
\author{Giovanni Petris, \email{GPetris@uark.edu}}

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

