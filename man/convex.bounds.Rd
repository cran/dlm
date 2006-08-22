\name{convex.bounds}
\alias{convex.bounds}
\title{Find the boundaries of a convex set}
\description{
Finds the boundaries of a bounded convex set along a specified 
straight line, using a bisection approach. It is mainly intended for 
use within \code{\link{arms}}.
}
\usage{
convex.bounds(x, dir, indFunc, ..., tol=1e-07)
}
\arguments{
  \item{x}{A point within the set}
  \item{dir}{A vector specifying a direction}
  \item{indFunc}{Indicator function of the set}
  \item{...}{Parameters passed to \code{indFunc}}
  \item{tol}{Tolerance}
}
\details{
Uses a bisection algorithm along a line having parametric representation
\code{x + t * dir}.
}
\value{
A vector \code{ans} of length two. The boundaries of the set are
\code{x + ans[1] * dir} and \code{x + ans[2] * dir}.
}
\author{Giovanni Petris <\email{GPetris@uark.edu}>}

\examples{
## boundaries of a unit circle
convex.bounds(c(0,0), c(1,1), indFunc=function(x) crossprod(x)<1)
}

\keyword{misc}% at least one, from doc/KEYWORDS

