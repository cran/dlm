\name{FF}
\alias{FF}
\alias{FF<-}
\alias{V}
\alias{V<-}
\alias{GG}
\alias{GG<-}
\alias{W}
\alias{W<-}
\alias{m0}
\alias{m0<-}
\alias{C0}
\alias{C0<-}

\title{Components of a dlm object}
\description{
  Functions to get or set specific components of an object of class \code{"dlm"}
}
\usage{
FF(x)
FF(x) <- value
V(x)
V(x) <- value
GG(x)
GG(x) <- value
W(x)
W(x) <- value
m0(x)
m0(x) <- value
C0(x)
C0(x) <- value
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{"dlm"}.}
  \item{value}{a numeric matrix (or vector for \code{m0}).} 
}
\details{
  Missing values are not allowed in \code{value}. The dimension of
  \code{value} must match the dimension of the current value of the
  specific component in \code{x}
}
\value{
  For the assignment forms, the updated \code{"dlm"} object.

  For the other forms, the specific component of \code{x}.
}
\author{Giovanni Petris \email{GPetris@uark.edu}}

\seealso{\code{\link{dlm}}}
\examples{
set.seed(222)
mod <- dlmRandom(5, 6)
all.equal( FF(mod), mod$FF )
all.equal( V(mod), mod$V )
all.equal( GG(mod), mod$GG )
all.equal( W(mod), mod$W )
all.equal( m0(mod), mod$m0 )
all.equal( C0(mod), mod$C0)
m0(mod)
m0(mod) <- rnorm(6)
C0
C0(mod) <- rwishart(10, 6)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}

