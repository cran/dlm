\name{dlmSmooth}
\alias{dlmSmooth}
\alias{dlmSmooth.dlmFiltered}
\alias{dlmSmooth.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{DLM smoothing}
\description{
  The function apply Kalman smoother to compute 
  smoothed values of the state vectors, together with their
  variance/covariance matrices. 
}
\usage{
dlmSmooth(y, \dots)
\method{dlmSmooth}{default}(y, mod, \dots)
\method{dlmSmooth}{dlmFiltered}(y, \dots, debug = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{an object used to select a method.}
  \item{\dots}{futher arguments passed to or from other methods.}
  \item{mod}{an object of class \code{"dlm"}.}
  \item{debug}{if \code{debug=FALSE}, faster C code will be used, otherwise
    all the computations will be performed in R.}
}
\details{%
  The default method returns means and variances of the smoothing
  distribution for a data vector (or matrix) \code{y} and a model
  \code{mod}.

  \code{dlmSmooth.dlmFiltered} produces the same output based on a
  \code{dlmFiltered} object, typically one produced by a call to
  \code{dlmFilter}. 

  The calculations are based on the singular value decomposition (SVD)
  of the relevant matrices. Variance matrices are returned in terms of
  their SVD.
}
\value{%
  A list with components
  \item{s}{Time series (or matrix) of smoothed values of the state
    vectors. The series starts one time unit before the first observation.} 
  \item{U.S}{See below.}
  \item{D.S}{Together with \code{U.S}, it gives the SVD of the variances
    of the smoothing errors.}
}
\section{Warning}{
  The observation variance \code{V} in \code{mod} must be nonsingular.
}
\references{Zhang, Y. and Li, X.R., Fixed-interval smoothing algorithm
  based on singular value decomposition, \emph{Proceedings of the 1996
    IEEE International Conference on Control Applications}.\cr
  Giovanni Petris (2010), An R Package for Dynamic Linear
  Models. Journal of Statistical Software, 36(12), 1-16.
  \url{https://www.jstatsoft.org/v36/i12/}.\cr
  Petris, Petrone, and Campagnoli, Dynamic Linear Models with
  R, Springer (2009).
}

\author{Giovanni Petris \email{GPetris@uark.edu}}

\seealso{See \code{\link{dlm}} for a description of dlm objects,
  \code{\link{dlmSvd2var}} to obtain a variance matrix from its SVD,
  \code{\link{dlmFilter}} for Kalman filtering, 
  \code{\link{dlmMLE}} for maximum likelihood estimation, and
  \code{\link{dlmBSample}} for drawing from the posterior distribution
  of the state vectors.} 

\examples{%
s <- dlmSmooth(Nile, dlmModPoly(1, dV = 15100, dW = 1470))
plot(Nile, type ='o')
lines(dropFirst(s$s), col = "red")

## Multivariate
set.seed(2)
tmp <- dlmRandom(3, 5, 20)
obs <- tmp$y
m <- tmp$mod
rm(tmp)

f <- dlmFilter(obs, m)
s <- dlmSmooth(f)
all.equal(s, dlmSmooth(obs, m))
}

\keyword{ts}
\keyword{smooth}
\keyword{misc}
