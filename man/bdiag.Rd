\name{bdiag}
\alias{bdiag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Build a block diagonal matrix}
\description{
  The function builds a block diagonal matrix. 
}
\usage{
bdiag(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{individual matrices, or a list of matrices.}
}
\details{
}
\value{
  A matrix obtained by combining the arguments.
}
\author{Giovanni Petris \email{GPetris@uark.edu}}

\examples{
bdiag(matrix(1:4,2,2),diag(3))
bdiag(matrix(1:6,3,2),matrix(11:16,2,3))
}
\keyword{misc}% at least one, from doc/KEYWORDS

