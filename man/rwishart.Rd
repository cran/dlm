\name{rwishart}
\alias{rwishart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Random Wishart matrix}
\description{
Generate a draw from a Wishart distribution.
}
\usage{
rwishart(df, p = nrow(SqrtSigma), Sigma, SqrtSigma = diag(p))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{degrees of freedom. It has to be integer.}
  \item{p}{dimension of the matrix to simulate.}
  \item{Sigma}{the matrix parameter Sigma of the Wishart distribution.}
  \item{SqrtSigma}{a \emph{square root} of the matrix parameter Sigma of the 
  Wishart distribution. Sigma must be equal to \code{crossprod(SqrtSigma)}.}
}
\details{
  The Wishart is a distribution on the set of nonnegative definite
  symmetric matrices. Its density is
  \deqn{p(W) = \frac{c |W|^{(n-p-1)/2}}{|\Sigma|^{n/2}}
    \exp\left\{-\frac{1}{2}\mathrm{tr}(\Sigma^{-1}W)\right\}}{%
    p(W) = c|W|^((n-p-1)/2) / |Sigma|^(n/2) exp(-tr(Sigma^(-1)W)/2)}
  where \eqn{n} is the degrees of freedom parameter \code{df} and
  \eqn{c} is a normalizing constant.
  The mean of the Wishart distribution is \eqn{n\Sigma}{n Sigma} and the
  variance of an entry is
  \deqn{\mathrm{Var}(W_{ij}) = n (\Sigma_{ij}^2 +
    \Sigma{ii}\Sigma_{jj})}{%
    Var(W[i,j]) = n (Sigma[i,j]^2 + Sigma[i,i] Sigma[j,j])}
  The matrix parameter, which should be a positive definite symmetric
  matrix, can be specified via either the argument Sigma or
  SqrtSigma. If Sigma is specified, then SqrtSigma is ignored. No checks
  are made for symmetry and positive definiteness of Sigma.
}
\value{
  The function returns one draw from the Wishart distribution with 
  \code{df} degrees of freedom and matrix parameter \code{Sigma} or
  \code{crossprod(SqrtSigma)}	
}
\references{Press (1982). Applied multivariate analysis.}
\author{Giovanni Petris \email{GPetris@uark.edu}}
\note{From a suggestion of B.Venables, posted on S-news}
\section{Warning}{The function only works for an integer number 
of degrees of freedom.}
 
\examples{
rwishart(25, p = 3)
a <- matrix(rnorm(9), 3)
rwishart(30, SqrtSigma = a)
b <- crossprod(a)
rwishart(30, Sigma = b)
}

\keyword{distribution}% at least one, from doc/KEYWORDS

