\name{dlmFilter}
\alias{dlmFilter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{DLM filtering}
\description{
  The functions applies Kalman filter to compute filtered
  values of the state vectors, together with their
  variance/covariance matrices. By default the function returns an object 
  of class \code{"dlmFiltered"}. Methods for \code{residuals} and \code{tsdiag} 
  for objects of class \code{"dlmFiltered"} exist.
}
\usage{
dlmFilter(y, mod, debug = FALSE, simplify = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{the data. \code{y} can be a vector, a matrix, a univariate or
    multivariate time series.}
  \item{mod}{an object of class \code{dlm}, or a list with components
    \code{m0}, \code{C0}, 
    \code{FF}, \code{V}, \code{GG}, \code{W}, and optionally \code{JFF},
    \code{JV}, \code{JGG}, \code{JW}, and \code{X}, defining the model
    and the parameters of the prior distribution.}
  \item{debug}{if \code{FALSE}, faster C code will be used, otherwise
    all the computations will be performed in R.}
  \item{simplify}{should the data be included in the output?}
}
\details{%
    The calculations are based on the singular value decomposition (SVD)
    of the relevant matrices. Variance matrices are returned in terms of
    their SVD.

    Missing values are allowed in \code{y}.
}
\value{%
  A list with the components described below. If \code{simplify} is 
  \code{FALSE}, the returned list has class \code{"dlmFiltered"}. 
  \item{y}{The input data, coerced to a matrix. This is present only if 
  \code{simplify} is \code{FALSE}.}
  \item{mod}{The argument \code{mod} (possibly simplified).}
  \item{m}{Time series (or matrix) of filtered values of the state vectors.
    The series starts one time unit before the first observation.}
  \item{U.C}{See below.}
  \item{D.C}{Together with \code{U.C}, it gives the SVD
    of the variances of the estimation errors.
    The variance of \eqn{m[t,]-theta[t,]}{m_t-\theta_t} is given by
    \code{U.C[[t]] \%*\% diag(D.C[t,]^2) \%*\% t(U.C[[t]])}.}
  \item{a}{Time series (or matrix) of predicted values of the state vectors
    given the observations up and including the previous time unit.}
  \item{U.R}{See below.}
  \item{D.R}{Together with \code{U.R}, it gives the SVD
    of the variances of the prediction errors.
    The variance of \eqn{a[t,]-theta[t,]}{a_t-\theta_t} is given by
    \code{U.R[[t]] \%*\% diag(D.R[t,]^2) \%*\% t(U.R[[t]])}.}
  \item{f}{Time series (or matrix) of one-step-ahead forecast of the
    observations.}
}
\section{Warning}{
  The observation variance \code{V} in \code{mod} must be nonsingular.
}
\references{Zhang, Y. and Li, X.R., Fixed-interval smoothing algorithm
  based on singular value decomposition, \emph{Proceedings of the 1996
    IEEE International Conference on Control Applications}.\cr
  Giovanni Petris (2010), An R Package for Dynamic Linear
  Models. Journal of Statistical Software, 36(12), 1-16.
  \url{http://www.jstatsoft.org/v36/i12/}.\cr
  Petris, Petrone, and Campagnoli, Dynamic Linear Models with R,
  Springer (2009).   
}
\author{Giovanni Petris \email{GPetris@uark.edu}}

\seealso{See \code{\link{dlm}} for a description of dlm objects,
  \code{\link{dlmSvd2var}} to obtain a variance matrix from its SVD,
  \code{\link{dlmMLE}} for maximum likelihood estimation, 
  \code{\link{dlmSmooth}} for Kalman smoothing, and
  \code{\link{dlmBSample}} for drawing from the posterior distribution
  of the state vectors.}

\examples{
nileBuild <- function(par) {
  dlmModPoly(1, dV = exp(par[1]), dW = exp(par[2]))
}
nileMLE <- dlmMLE(Nile, rep(0,2), nileBuild); nileMLE$conv
nileMod <- nileBuild(nileMLE$par)
V(nileMod)
W(nileMod)
nileFilt <- dlmFilter(Nile, nileMod)
nileSmooth <- dlmSmooth(nileFilt)
plot(cbind(Nile, nileFilt$m[-1], nileSmooth$s[-1]), plot.type='s',
     col=c("black","red","blue"), ylab="Level", main="Nile river", lwd=c(1,2,2))
}

\keyword{ts}
\keyword{misc}
