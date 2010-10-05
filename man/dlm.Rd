\name{dlm}
\alias{dlm}
\alias{as.dlm}
\alias{is.dlm}
\title{dlm objects}
\description{
  The function \code{dlm} is used to create Dynamic Linear Model objects.
  \code{as.dlm} and \code{is.dlm} coerce an object to a Dynamic Linear
  Model object and test whether an object is a Dynamic Linear Model. 
}
\usage{
dlm(...)
as.dlm(obj)
is.dlm(obj)
}
\arguments{
  \item{...}{list with named elements \code{m0},
    \code{C0}, \code{FF}, \code{V}, \code{GG}, \code{W} and,
    optionally, 
    \code{JFF}, \code{JV}, \code{JGG}, \code{JW}, and
    \code{X}. The first six are the usual vector and matrices that define
    a time-invariant DLM. The remaining elements are used for time-varying
    DLM. \code{X}, if present, should be a matrix. If \code{JFF} is not
    \code{NULL}, then it must be a matrix of the 
    same dimension of \code{FF}, with the \eqn{(i,j)} element being zero if
    \code{FF[i,j]} is time-invariant, and a positive integer \eqn{k}
    otherwise. In this case the \eqn{(i,j)} element of \code{FF} at time
    \eqn{t} will be \code{X[t,k]}. A similar interpretation holds for
    \code{JV}, \code{JGG}, and \code{JW}. \code{...} may have additional
    components, that are not used by \code{dlm}. The named components
    may also be passed to the function as individual arguments.}     
  \item{obj}{an arbitrary \R object.}
}
\details{
  The function \code{dlm} is used to create Dynamic Linear Model
  objects. These are lists with the named elements described above and
  with class of \code{"dlm"}. 

  Class \code{"dlm"} has a number of methods. In particular, consistent
  DLM can be added together to produce another DLM.
}
\value{
  For \code{dlm}, an object of class \code{"dlm"}.}
\references{Giovanni Petris (2010), An R Package for Dynamic Linear
  Models. Journal of Statistical Software, 36(12), 1-16.
  \url{http://www.jstatsoft.org/v36/i12/}.\cr
  Petris, Petrone, and Campagnoli, Dynamic Linear Models with
  R, Springer (2009).\cr
  West and Harrison, Bayesian forecasting and
  dynamic models (2nd ed.), Springer (1997).
} 
\seealso{\code{\link{dlmModReg}}, \code{\link{dlmModPoly}},
  \code{\link{dlmModARMA}}, \code{\link{dlmModSeas}}, to create
  particular objects of class \code{"dlm"}.   
}
\examples{
## Linear regression as a DLM
x <- matrix(rnorm(10),nc=2)
mod <- dlmModReg(x)
is.dlm(mod)

## Adding dlm's
dlmModPoly() + dlmModSeas(4) # linear trend plus quarterly seasonal component
}
\author{Giovanni Petris \email{GPetris@uark.edu}}
\keyword{misc}
