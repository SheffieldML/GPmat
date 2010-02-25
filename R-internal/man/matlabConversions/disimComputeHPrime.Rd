\name{disimComputeHPrime}
\Rdversion{1.0}
\alias{disimComputeHPrime}
\title{Helper function for comptuing part of the DISIM kernel.}
\description{
  computes a portion of the SIM kernel and gradients with respect to various parameters.
}
\usage{
values <- disimComputeHPrime(t1, t2, decay0, decay1, decay2, l,
                             t1, t2, decay0, decay1, decay2, l)
}
\arguments{
  \item{t1}{first time input (number of time points x 1).}
  \item{t2}{second time input (number of time points x 1).}
  \item{decay0}{Decay rate for the driving system.}
  \item{decay1}{Decay rate for first system.}
  \item{decay2}{Decay rate for second system.}
  \item{l}{length scale of latent process.}
  \item{t1}{first time input (number of time points x 1).}
  \item{t2}{second time input (number of time points x 1).}
  \item{decay0}{Decay rate for the driving system.}
  \item{decay1}{Decay rate for first system.}
  \item{decay2}{Decay rate for second system.}
  \item{l}{length scale of latent process.}
}
\value{
  \item{h}{result of this subcomponent of the kernel for the given values.}
  \item{h}{result of this subcomponent of the kernel for the given values.}
  \item{grad_D_decay0}{gradient of H with respect to DECAY0.}
  \item{grad_D_decay1}{gradient of H with respect to DECAY1.}
  \item{grad_D_decay2}{gradient of H with respect to DECAY2.}
  \item{grad_L}{gradient of H with respect to length scale of latent process.}
}
\seealso{
\code{\link{disimKernParamInit}, \link{disimComputeH}}.
}
\examples{
## missing
}
\keyword{model}
