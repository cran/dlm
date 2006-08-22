\name{rwishart}
\alias{rwishart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Random Whishart matrix}
\description{
Generate a draw from a Wishart distribution.
}
\usage{
rwishart(df, p = nrow(SqrtSigma), SqrtSigma = diag(p))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{Degrees of freedom. It has to be integer.}
  \item{p}{Dimension of the matrix to simulate.}
  \item{SqrtSigma}{A "square root" of the matrix parameter Sigma of the 
  Wishart distribution. Sigma must be equal to \code{crossprod(SqrtSigma)}.}
}
\value{
  The function returns one draw from the Wishart distribution with 
  \code{df} degrees of freedom and matrix parameter 
  \code{crossprod(SqrtSigma)}	
}
\references{Mardia, Kent and Bibby (1979). Multivariate analysis.}
\author{B.Venables, from S-news.}
\note{}
\section{Warning}{The function only works for an integer number 
of degrees of freedom.}
 
\examples{
rwishart(25, p = 2)
}

\keyword{distribution}% at least one, from doc/KEYWORDS

