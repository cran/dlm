\name{dlmBSample}
\alias{dlmBSample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Draw from the posterior distribution of the state vectors}
\description{
  The function simulates one draw from the posterior distribution
  of the state vectors. 
}
\usage{
dlmBSample(modFilt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{modFilt}{%
    a list, typically the ouptut from \code{dlmFilter}, with elements
    \code{m}, \code{U.C}, \code{D.C}, \code{a}, \code{U.R}, \code{D.R} (see
    the value returned by \code{dlmFilter}), and \code{mod} The latter is an
    object of class \code{"dlm"} or a list with elements \code{GG}, \code{W}
    and, optionally, \code{JGG}, \code{JW}, and \code{X}}
}
\details{
  The calculations are based on singular value decomposition.
}
\value{
  The function returns a draw from the posterior distribution 
  of the state vectors. If \code{m} is a time series then the returned 
  value is a time series with the same \code{tsp}, otherwise it is 
  a matrix or vector.
}
\references{Petris, Petrone, and Campagnoli, Dynamic Linear Models with
  R, Springer (2009). West and Harrison, Bayesian forecasting and
  dynamic models (2nd ed.), Springer (1997).
} 
\author{Giovanni Petris \email{GPetris@uark.edu}}

\seealso{See also \code{\link{dlmFilter}}}
\examples{
nileMod <- dlmModPoly(1, dV = 15099.8, dW = 1468.4)
nileFilt <- dlmFilter(Nile, nileMod)
nileSmooth <- dlmSmooth(nileFilt) # estimated "true" level
plot(cbind(Nile, nileSmooth$s[-1]), plot.type = "s",
     col = c("black", "red"), ylab = "Level",
     main = "Nile river", lwd = c(2, 2)) 
for (i in 1:10) # 10 simulated "true" levels 
    lines(dlmBSample(nileFilt[-1]), lty=2) 
}
\keyword{misc}% at least one, from doc/KEYWORDS

