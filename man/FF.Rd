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
\alias{FF.dlm}
\alias{FF<-.dlm}
\alias{V.dlm}
\alias{V<-.dlm}
\alias{GG.dlm}
\alias{GG<-.dlm}
\alias{W.dlm}
\alias{W<-.dlm}
\alias{m0.dlm}
\alias{m0<-.dlm}
\alias{C0.dlm}
\alias{C0<-.dlm}
\alias{JFF}
\alias{JFF<-}
\alias{JV}
\alias{JV<-}
\alias{JGG}
\alias{JGG<-}
\alias{JW}
\alias{JW<-}
\alias{X}
\alias{X<-}
\alias{JFF.dlm}
\alias{JFF<-.dlm}
\alias{JV.dlm}
\alias{JV<-.dlm}
\alias{JGG.dlm}
\alias{JGG<-.dlm}
\alias{JW.dlm}
\alias{JW<-.dlm}
\alias{X.dlm}
\alias{X<-.dlm}


\title{Components of a dlm object}
\description{
  Functions to get or set specific components of an object of class \code{dlm}
}
\usage{
\method{FF}{dlm}(x)
\method{FF}{dlm}(x) <- value
\method{V}{dlm}(x)
\method{V}{dlm}(x) <- value
\method{GG}{dlm}(x)
\method{GG}{dlm}(x) <- value
\method{W}{dlm}(x)
\method{W}{dlm}(x) <- value
\method{m0}{dlm}(x)
\method{m0}{dlm}(x) <- value
\method{C0}{dlm}(x)
\method{C0}{dlm}(x) <- value
\method{JFF}{dlm}(x)
\method{JFF}{dlm}(x) <- value
\method{JV}{dlm}(x)
\method{JV}{dlm}(x) <- value
\method{JGG}{dlm}(x)
\method{JGG}{dlm}(x) <- value
\method{JW}{dlm}(x)
\method{JW}{dlm}(x) <- value
\method{X}{dlm}(x)
\method{X}{dlm}(x) <- value
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{dlm}.}
  \item{value}{a numeric matrix (or vector for \code{m0}).} 
}
\details{
  Missing or infinite values are not allowed in \code{value}. The dimension of
  \code{value} must match the dimension of the current value of the
  specific component in \code{x}
}
\value{
  For the assignment forms, the updated \code{dlm} object.

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
C0(mod)
C0(mod) <- rwishart(10, 6)
### A time-varying model
mod <- dlmModReg(matrix(rnorm(10), 5, 2))
JFF(mod)
X(mod)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}

