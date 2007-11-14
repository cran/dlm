\name{dlmSum}
\alias{dlmSum}
\alias{\%+\%}
\title{Outer sum of Dynamic Linear Models}
\description{
  \code{dlmSum} creates a unique DLM out of two or more
  independent DLMs. \code{\%+\%} is an alias for \code{dlmSum}.
}
\usage{
dlmSum(...)
x \%+\% y
}

\arguments{
  \item{\dots}{any number of objects of class \code{dlm}, or a list of
    such objects.}
  \item{x, y}{objects of class \code{dlm}.}
}
\value{
  An object of class \code{dlm}, representing the outer sum of the
  arguments.
}

\author{Giovanni Petris \email{GPetris@uark.edu}}

\examples{
m1 <- dlmModPoly(2)
m2 <- dlmModPoly(1)
dlmSum(m1, m2)
m1 \%+\% m2 # same thing
}

\keyword{ts}
\keyword{misc}
