\name{dlmGibbsDIG}
\alias{dlmGibbsDIG}
\title{Gibbs sampling for d-inverse-gamma model}
\description{
  The function implements a Gibbs sampler for a univariate DLM
  having one or more unknown variances in its specification.
}
\usage{
dlmGibbsDIG(y, mod, a.y, b.y, a.theta, b.theta, shape.y, rate.y,
            shape.theta, rate.theta, n.sample = 1,
            thin = 0, ind, save.states = TRUE,
            progressBar = interactive())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{data vector or univariate time series}
  \item{mod}{a dlm for univariate observations}
  \item{a.y}{prior mean of observation variance}
  \item{b.y}{prior variance of observation variance}
  \item{a.theta}{prior mean of system variances}
  \item{b.theta}{prior variance of system variances}
  \item{shape.y}{shape parameter of the prior of observation variance}
  \item{rate.y}{rate parameter of the prior of observation variance}
  \item{shape.theta}{shape parameter of the prior of system variances}
  \item{rate.theta}{rate parameter of the prior of system variances}
  \item{n.sample}{requested number of Gibbs iterations}
  \item{thin}{discard \code{thin} iterations for every saved iteration}
  \item{ind}{indicator of the system variances that need to be estimated}
  \item{save.states}{should the simulated states be included in the output?}
  \item{progressBar}{should a text progress bar be displayed during execution?}
}
\details{
  The \emph{d-inverse-gamma} model is a constant univariate DLM with unknown
  observation variance, diagonal system variance with unknown diagonal
  entries.  Some of these entries may be known, in which case they are
  typically zero.  Independent inverse gamma priors are assumed for the
  unknown variances.  These can be specified be mean and variance or,
  alternatively, by shape and rate.  Recycling is applied for the prior
  parameters of unknown system variances.  The argument \code{ind} can
  be used to specify the index of the unknown system variances, in case
  some of the diagonal elements of \code{W} are known.  The unobservable
  states are generated in the Gibbs sampler and are returned if
  \code{save.states = TRUE}.  For more details on the model and usage
  examples, see the package vignette. 
}
\value{
  The function returns a list of simulated values.
  \item{dV}{simulated values of the observation variance.}
  \item{dW}{simulated values of the unknown diagonal elements of the
    system variance.}
  \item{theta}{simulated values of the state vectors.}
}

\author{Giovanni Petris \email{GPetris@uark.edu}}

\examples{
## See the package vignette for an example
}

\keyword{misc}

