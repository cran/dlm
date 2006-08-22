\name{dlmSmooth}
\alias{dlmSmooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{DLM smoothing}
\description{
  The function apply Kalman smoother to compute 
  smoothed values of the state vectors, together with their
  variance/covariance matrices. 
}
\usage{
dlmSmooth(modFilt, debug = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modFilt}{An object of class \code{"dlmFiltered"}, such as returned by 
    \code{dlmFilter}}
  \item{debug}{If \code{FALSE}, faster C code will be used, otherwise
    all the computations will be performed in R.}
}
\details{%
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
\references{Zhang, Y. and Li, X.R., Fixed-interval smoothing algorithm
  based on singular value decomposition, \emph{Proceedings of the 1996
    IEEE International Conference on Control Applications}.}
\author{Giovanni Petris, \email{<GPetris@uark.edu>}}

\seealso{See \code{\link{dlm}} for a description of dlm objects,
  \code{\link{dlmSvd2var}} to obtain a variance matrix from its SVD,
  \code{\link{dlmFilter}} for Kalman filtering, 
  \code{\link{dlmMLE}} for maximum likelihood estimation, and
  \code{\link{dlmBSample}} for drawing from the posterior distribution
  of the state vectors.} 

\examples{%
### See example for dlmFilter
}

\keyword{ts}
\keyword{smooth}
\keyword{misc}
