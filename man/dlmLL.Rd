\name{dlmLL}
\alias{dlmLL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Log likelihood evaluation for a state space model}
\description{
	Function that computes the log likelihood of a state space model. 
}
\usage{
dlmLL(y, mod, debug=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{a vector, matrix, or time series of data.}
  \item{mod}{an object of class \code{"dlm"}, or a list with components
    \code{m0}, \code{C0}, 
    \code{FF}, \code{V}, \code{GG}, \code{W} defining the model
    and the parameters of the prior distribution.} 
  \item{debug}{if \code{debug=TRUE}, the function uses R code, otherwise
    it uses faster C code.}
}
\details{
  The calculations are based on singular value decomposition.
  Missing values are allowed in \code{y}.
}
\value{
  The function returns the negative of the loglikelihood.
}
\section{Warning}{
  The observation variance \code{V} in \code{mod} must be nonsingular.
}
\references{Durbin and Koopman, Time series analysis by state space
  methods, Oxford University Press, 2001.}
\author{Giovanni Petris \email{GPetris@uark.edu}}


\seealso{\code{\link{dlmMLE}}, \code{\link{dlmFilter}} for the definition of 
  the equations of the model.}
\examples{
##---- See the examples for dlmMLE ----
}
\keyword{misc}% at least one, from doc/KEYWORDS

