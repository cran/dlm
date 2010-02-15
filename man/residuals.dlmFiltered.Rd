\name{residuals.dlmFiltered}
\alias{residuals.dlmFiltered}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{One-step forecast errors}
\description{
  The function computes one-step forecast errors for a filtered dynamic
  linear model.
}
\usage{
\method{residuals}{dlmFiltered}(object, ..., type = c("standardized", "raw"), sd = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{"dlmFiltered"}, such as the
    output from \code{dlmFilter}}
  \item{\dots}{unused additional arguments.}
  \item{type}{should standardized or raw forecast errors be produced?}
  \item{sd}{when \code{sd = TRUE}, standard deviations are returned as well.}
}
\value{
  A vector or matrix (in the multivariate case) of one-step forecast
  errors, standardized if \code{type = "standardized"}. Time series
  attributes of the original observation vector (matrix) are retained by
  the one-step forecast errors.

  If \code{sd = TRUE} then the returned value is a list with the
  one-step forecast errors in component \code{res} and the corresponding
  standard deviations in component \code{sd}. 
}
\note{The \code{object} argument must include a component \code{y}
  containing the data. This component will not be present if
  \code{object} was obtained by calling \code{dlmFilter} with
  \code{simplify = TRUE}.}
\references{Petris, Petrone, and Campagnoli, Dynamic Linear Models with
  R, Springer (2009). Harrison and West, Bayesian forecasting and dynamic models
  (2nd ed.), Springer (1997).}
\author{Giovanni Petris \email{GPetris@uark.edu}}

\seealso{\code{\link{dlmFilter}}}
\examples{
## diagnostic plots 
nileMod <- dlmModPoly(1, dV = 15100, dW = 1468)
nileFilt <- dlmFilter(Nile, nileMod)
res <- residuals(nileFilt, sd=FALSE)
qqnorm(res)
tsdiag(nileFilt)
}
\keyword{misc}% at least one, from doc/KEYWORDS

