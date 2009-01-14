\name{dlmModReg}
\alias{dlmModReg}
\title{Create a DLM representation of a regression model}
\description{
  The function creates a dlm representation of a linear regression
  model. 
}
\usage{
dlmModReg(X, addInt = TRUE, dV = 1, dW = rep(0, NCOL(X) + addInt),
          m0 = rep(0, length(dW)),
          C0 = 1e+07 * diag(nrow = length(dW)))
}
\arguments{
  \item{X}{the design matrix}
  \item{addInt}{logical: should an intercept be added?}
  \item{dV}{variance of the observation noise.}
  \item{dW}{diagonal elements of the variance matrix of the system noise.}
  \item{m0}{\eqn{m_0}{m0}, the expected value of the pre-sample state vector.}
  \item{C0}{\eqn{C_0}{C0}, the variance matrix of the pre-sample state vector.}
}
\details{
  By setting \code{dW} equal to a nonzero vector one obtains a DLM
  representation of a dynamic regression model. The default value zero
  of \code{dW} corresponds to standard linear regression. Only
  univariate regression is currently covered. 
}
\value{
  An object of class dlm representing the specified regression model.
}
\references{West and Harrison, Bayesian forecasting and dynamic models
  (2nd ed.), Springer, 1997.}
\author{Giovanni Petris, \email{GPetris@uark.edu}}

\seealso{\code{\link{dlmModARMA}}, \code{\link{dlmModPoly}},
  \code{\link{dlmModSeas}}}
\examples{
x <- matrix(runif(6,4,10), nc = 2); x
dlmModReg(x)
dlmModReg(x, addInt = FALSE)
}
\keyword{misc}% at least one, from doc/KEYWORDS

