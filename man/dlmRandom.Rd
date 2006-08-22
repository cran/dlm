\name{dlmRandom}
\alias{dlmRandom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Random DLM}
\description{
  Generate a random (constant or time-varying) object of class
  \code{"dlm"}, together with states and observations from it.
}
\usage{
dlmRandom(m, p, nobs = 0, JFF, JV, JGG, JW)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{Dimension of the observation vector.}
  \item{p}{Dimension of the state vector.}
  \item{nobs}{Number of states and observations to simulate from the model.}
  \item{JFF}{Should the model have a time-varying \code{FF} component?}
  \item{JV}{Should the model have a time-varying \code{V} component?}
  \item{JGG}{Should the model have a time-varying \code{GG} component?}
  \item{JW}{Should the model have a time-varying \code{W} component?}
}
\details{
  The function generates randomly the system and oservation matrices and
  the variances of a DLM having the specified state and observation
  dimension. The system matrix \code{GG} is guaranteed to have
  eigenvalues strictly less than one, which means that a constant DLM is
  asymptotically stationary. The default behavior is to generate a
  constant DLM. If \code{JFF} is \code{TRUE} then a model for
  \code{nobs} observations in which all
  the elements of \code{FF} are time-varying is generated. Similarly
  with \code{JV}, \code{JGG}, and \code{JW}.
}
\value{
  The function returns a list with the following components.
  \item{mod}{An object of class \code{"dlm"}.}
  \item{theta}{Matrix of simulated state vectors from the model.}
  \item{y}{Matrix of simulated observations from the model.}
  If \code{nobs} is zero, only the \code{mod} component is returned.
}
\references{Anderson and Moore, Optimal filtering, Prentice-Hall (1979)}
\author{Giovanni Petris \email{<GPetris@uark.edu>}}
\seealso{\code{\link{dlm}}}
\examples{
dlmRandom(1, 3, 5)
}
\keyword{misc}% at least one, from doc/KEYWORDS
\keyword{datagen}% __ONLY ONE__ keyword per line
