\name{dlmModARMA}
\alias{dlmModARMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create a DLM representation of an ARMA process}
\description{
  The function creates an object of class dlm representing a specified
  univariate or multivariate ARMA process
}
\usage{
dlmModARMA(ar = NULL, ma = NULL, sigma2 = 1, dV, m0, C0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ar}{a vector or a list of matrices (in the multivariate case)
    containing the autoregressive coefficients.}
  \item{ma}{a vector or a list of matrices (in the multivariate case)
    containing the moving average coefficients.}
  \item{sigma2}{the variance (or variance matrix) of the innovations.}
  \item{dV}{the variance, or the diagonal elements of the variance
    matrix in the multivariate case, of the observation noise. \code{V}
    is assumed to be diagonal and it defaults to zero.}
  \item{m0}{\eqn{m_0}{m0}, the expected value of the pre-sample state vector.}
  \item{C0}{\eqn{C_0}{C0}, the variance matrix of the pre-sample state vector.}
}
\details{
  The returned DLM only gives one of the many possible representations
  of an ARMA process. 
%  What is the name of this representation???
}
\value{
  The function returns an object of class dlm representing the ARMA
  model specified by \code{ar}, \code{ma}, and \code{sigma2}.
}
\references{Giovanni Petris (2010), An R Package for Dynamic Linear
  Models. Journal of Statistical Software, 36(12), 1-16.
  \url{http://www.jstatsoft.org/v36/i12/}.\cr
  Petris, Petrone, and Campagnoli, Dynamic Linear Models with
  R, Springer (2009).\cr
  Durbin and Koopman, Time series analysis by state space 
  methods, Oxford University Press, 2001.
}

\author{Giovanni Petris \email{GPetris@uark.edu}}

\seealso{\code{\link{dlmModPoly}}, \code{\link{dlmModSeas}},
  \code{\link{dlmModReg}}}
\examples{
## ARMA(2,3)
dlmModARMA(ar = c(.5,.1), ma = c(.4,2,.3), sigma2=1)
## Bivariate ARMA(2,1)
dlmModARMA(ar = list(matrix(1:4,2,2), matrix(101:104,2,2)),
           ma = list(matrix(-4:-1,2,2)), sigma2 = diag(2))
}
\keyword{misc}% at least one, from doc/KEYWORDS

