\name{dlmSvd2var}
\alias{dlmSvd2var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute a nonnegative definite matrix from its
  Singular Value Decomposition} 
\description{
  The function computes a nonnegative definite matrix from its Singular
  Value Decomposition. 
}
\usage{
dlmSvd2var(u, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{A square matrix, or a list of square matrices for a
    vectorized usage.}
  \item{d}{A vector, or a matrix for a vectorized usage.}
}
\details{
  The SVD of a nonnegative definite \eqn{n} by \eqn{n} square matrix
  \eqn{x} can be written as \eqn{u d^2 u'}, where \eqn{u} is an \eqn{n}
  by \eqn{n} orthogonal matrix and \eqn{d} is a diagonal matrix. For a
  single matrix, the function returns just \eqn{u d^2 u'}. Note that the
  argument \code{d} is a vector containing the diagonal elements of
  \eqn{d}. For a vectorized usage, \code{u} is a list of square
  matrices, and \code{d} is a matrix. The returned value in this case is
  a list of matrices, with the element \eqn{i} being \code{u[[i]] \%*\%
    diag(d[i,]^2) \%*\% t(u[[i]])}. 
}
\value{
  The function returns a nonnegative definite matrix, reconstructed from
  its SVD, or a list of such matrices (see details above).
}
\references{Horn and Johnson, Matrix analysis, Cambridge University
  Press (1985)}
\author{Giovanni Petris \email{GPetris@uark.edu}}

\examples{
x <- matrix(rnorm(16),4,4)
x <- crossprod(x)
tmp <- La.svd(x)
all.equal(dlmSvd2var(tmp$u, sqrt(tmp$d)), x)
## Vectorized usage
x <- dlmFilter(Nile, dlmModPoly(1, dV=15099, dW=1469))
x$se <- sqrt(unlist(dlmSvd2var(x$U.C, x$D.C)))
## Level with 50\% probability interval
plot(Nile, lty=2)
lines(window(x$m, start=start(Nile)), col="blue")
lines(window(x$m - .67*x$se, start=start(Nile)), lty=3, col="blue")
lines(window(x$m + .67*x$se, start=start(Nile)), lty=3, col="blue")
}
\keyword{array}% at least one, from doc/KEYWORDS
\keyword{misc}% __ONLY ONE__ keyword per line
