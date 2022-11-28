\name{mcmc}
\alias{mcmcMean}
\alias{mcmcMeans}
\alias{mcmcSD}
\alias{ergMean}

\title{Utility functions for MCMC output analysis}
\description{
  Returns the mean, the standard deviation of the mean, and a sequence
  of partial means of the input vector or matrix. 
}
\usage{
mcmcMean(x, sd = TRUE)
mcmcMeans(x, sd = TRUE)
mcmcSD(x)
ergMean(x, m = 1)
}

\arguments{
  \item{x}{vector or matrix containing the output of a Markov chain
    Monte Carlo simulation.}
  \item{sd}{logical: should an estimate of the Monte Carlo standard
    deviation be reported?} 
  \item{m}{ergodic means are computed for \code{i} in \code{m:NROW(x)}}
}

\details{
  The argument \code{x} is typically the output from a simulation. If a
  matrix, rows are considered consecutive simulations of a target
  vector. In this case means, standard deviations, and ergodic means
  are returned for each column. The standard deviation of the mean is
  estimated using Sokal's method (see the reference). \code{mcmcMeans}
  is an alias for \code{mcmcMean}.
}

\value{
  \code{mcmcMean} returns the sample mean of a vector containing the output
  of an MCMC sampler, together with an estimated standard error. If the input
  is a matrix, means and standard errors are computed for each column.\cr
  \code{mcmcSD} returns an estimate of the standard deviation of the mean for
  the output of an MCMC sampler.\cr
  \code{ergMean} returns a vector of running ergodic means.
}
  
\references{P. Green (2001). A Primer on Markov Chain Monte Carlo. In
  \emph{Complex Stochastic Systems}, (Barndorff-Nielsen, Cox and
  Kl\"uppelberg, eds.). Chapman and Hall/CRC. }
\author{Giovanni Petris \email{GPetris@uark.edu}}

\examples{
x <- matrix(rexp(1000), nc=4)
dimnames(x) <- list(NULL, LETTERS[1:NCOL(x)])
mcmcSD(x)
mcmcMean(x)
em <- ergMean(x, m = 51)
plot(ts(em, start=51), xlab="Iteration", main="Ergodic means")
}

\keyword{misc}

