\name{dlmModARMA}
\alias{dlmModARMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create a DLM representation of an ARMA process}
\description{
  The function creates an object of class dlm representing a specified
  univariate or multivariate ARMA process
}
\usage{
dlmModARMA(ar = NULL, ma = NULL, sigma2 = 1, m0, C0, dV)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ar}{a vector or a list of matrices (in the multivariate case)
    containing the autoregressive coefficients.}
  \item{ma}{a vector or a list of matrices (in the multivariate case)
    containing the moving average coefficients.}
  \item{sigma2}{the variance (or variance matrix) of the innovations.}
  \item{m0}{returned as \code{m0} in the created dlm object. It defaults
    to a vector of zeroes.}
  \item{C0}{returned as \code{C0} in the created dlm object. It defaults
    to a \code{1e7 * diag(k)}, where \code{k} is the dimension of the
    state vector in the DLM representation of the specified ARMA model.}
  \item{dV}{the variance, or the diagonal elements of the variance
    matrix, of the observational noise. \code{V} is assumed to be
    diagonal and it defaults to zero.}
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
\references{Durbin and Koopman, Time series analysis by state space 
  methods, Oxford University Press, 2001.}
\author{Giovanni Petris, \email{GPetris@uark.edu}}

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

