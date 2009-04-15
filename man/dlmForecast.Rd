\name{dlmForecast}
\alias{dlmForecast}
\title{Prediction and simulation of future observations}
\description{
  The function evaluates the expected value and variance of future
  observations and system states. It can also generate a sample from the
  distribution of future observations and system states.
}
\usage{
  dlmForecast(mod, nAhead = 1, method = c("plain", "svd"), sampleNew = FALSE)
}
\arguments{
  \item{mod}{an object of class \code{"dlm"}, or a list with components
    \code{m0}, \code{C0},
    \code{FF}, \code{V}, \code{GG}, and \code{W}, defining the model
    and the parameters of the prior distribution. \code{mod} can also be
    an object of class \code{"dlmFiltered"}, such as the output from
    \code{dlmFilter}.}
  \item{nAhead}{number of steps ahead for which a forecast is
    requested.}
  \item{method}{\code{method="svd"} uses singular value decomposition
    for the calculations. Currently, only \code{method="plain"}
    is implemented.}
  \item{sampleNew}{if \code{sampleNew=n} for an integer \code{n}, them a
    sample of size \code{n} from the forecast distribution of states and
    observables will be returned.}
}
\value{
  A list with components
  \tabular{ll}{
    \code{a}\tab matrix of expected values of future states\cr
    \code{R}\tab list of variances of future states\cr
    \code{f}\tab matrix of expected values of future observations\cr
    \code{Q}\tab list of variances of future observations\cr
    \code{newStates}\tab list of matrices containing the simulated
    future values\cr
      \tab of the states. Each component of the list corresponds\cr
      \tab to one simulation.\cr
    \code{newObs}\tab same as \code{newStates}, but for the observations.
  }
The last two components are not present if \code{sampleNew=FALSE}.
}
\note{
  The function is currently entirely written in \R and is not
  particularly fast. Currently, only constant models are allowed. 
}
\examples{
## Comparing theoretical prediction intervals with sample quantiles
set.seed(353)
n <- 20; m <- 1; p <- 5
mod <- dlmModPoly() + dlmModSeas(4, dV=0)
W(mod) <- rwishart(2*p,p) * 1e-1
m0(mod) <- rnorm(p, sd=5)
C0(mod) <- diag(p) * 1e-1
new <- 100
fore <- dlmForecast(mod, nAhead=n, sampleNew=new)
ciTheory <- (outer(sapply(fore$Q, FUN=function(x) sqrt(diag(x))), qnorm(c(0.1,0.9))) +
             as.vector(t(fore$f)))
ciSample <- t(apply(array(unlist(fore$newObs), dim=c(n,m,new))[,1,], 1,
                    FUN=function(x) quantile(x, c(0.1,0.9))))
plot.ts(cbind(ciTheory,fore$f[,1]),plot.type="s", col=c("red","red","green"),ylab="y")
for (j in 1:2) lines(ciSample[,j], col="blue")
legend(2,-40,legend=c("forecast mean", "theoretical bounds", "Monte Carlo bounds"),
       col=c("green","red","blue"), lty=1, bty="n")
}
\author{Giovanni Petris \email{GPetris@uark.edu}}
\keyword{misc}
