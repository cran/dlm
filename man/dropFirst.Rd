\name{dropFirst}
\alias{dropFirst}
\title{Drop the first element of a vector or matrix}
\description{
  A utility function, \code{dropFirst} drops the first element of a
  vector or matrix, retaining the correct time series attributes, in
  case the argument is a time series object. 
}
\usage{
dropFirst(x)
}
\arguments{
  \item{x}{a vector or matrix.}
}
\value{
  The function returns \code{x[-1]} or \code{x[-1,]}, if the argument is
  a matrix. For an argument of class \code{ts} the class is preserved,
  together with the correct \code{tsp} attribute.
}
\author{Giovanni Petris \email{GPetris@uark.edu}}

\examples{
(pres <- dropFirst(presidents))
start(presidents)
start(pres)
}


\keyword{misc}
\keyword{ts}
