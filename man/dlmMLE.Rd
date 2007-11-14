\name{dlmMLE}
\alias{dlmMLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Parameter estimation by maximum likelihood}
\description{
  The function returns the MLE of unknown parameters in the specification of 
  a state space model. 
}
\usage{
dlmMLE(y, parm, build, method = "L-BFGS-B", ..., debug = FALSE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{a vector, matrix, or time series of data.}
  \item{parm}{vector of initial values - for the optimization routine -
    of the unknown parameters.} 
  \item{build}{a function that takes a vector of the same length as
    \code{parm} and returns an object of class dlm, or a list that may
    be interpreted as such.}
  \item{method}{passed to \code{optim}.}
  \item{...}{additional arguments passed to \code{optim} and
    \code{build}.}
  \item{debug}{if \code{debug=TRUE}, the likelihood calculations are done
    entirely in R, otherwise C functions are used.}
}
\details{
  The evaluation of the loglikelihood is done by \code{dlmLL}. 
  For the optimization, \code{optim} is called. It is possible for the
  model to depend on additional parameters, other than those in
  \code{parm}, passed to \code{build} via the \code{...} argument.
}
\value{
  The function \code{dlmMLE} returns the value returned by \code{optim}.
}
\references{}
\author{Giovanni Petris \email{GPetris@uark.edu}}

\seealso{
  \code{\link{dlmLL}}, \code{\link{dlm}}.
}
  
\examples{
data(NelPlo)
### multivariate local level -- seemingly unrelated time series
buildSu <- function(x) {
  Vsd <- exp(x[1:2])
  Vcorr <- tanh(x[3])
  V <- Vsd \%o\% Vsd
  V[1,2] <- V[2,1] <- V[1,2] * Vcorr
  Wsd <- exp(x[4:5])
  Wcorr <- tanh(x[6])
  W <- Wsd \%o\% Wsd
  W[1,2] <- W[2,1] <- W[1,2] * Wcorr
  return(list(
              m0 = rep(0,2),
              C0 = 1e7 * diag(2),
              FF = diag(2),
              GG = diag(2),
              V = V,
              W = W))
}

suMLE <- dlmMLE(NelPlo, rep(0,6), buildSu); suMLE
buildSu(suMLE$par)[c("V","W")]
StructTS(NelPlo[,1], type="level") ## compare with W[1,1] and V[1,1]
StructTS(NelPlo[,2], type="level") ## compare with W[2,2] and V[2,2]

## multivariate local level model with homogeneity restriction
buildHo <- function(x) {
  Vsd <- exp(x[1:2])
  Vcorr <- tanh(x[3])
  V <- Vsd \%o\% Vsd
  V[1,2] <- V[2,1] <- V[1,2] * Vcorr
  return(list(
              m0 = rep(0,2),
              C0 = 1e7 * diag(2),
              FF = diag(2),
              GG = diag(2),
              V = V,
              W = x[4]^2 * V))
}

hoMLE <- dlmMLE(NelPlo, rep(0,4), buildHo); hoMLE
buildHo(hoMLE$par)[c("V","W")]
}
\keyword{misc}% at least one, from doc/KEYWORDS

