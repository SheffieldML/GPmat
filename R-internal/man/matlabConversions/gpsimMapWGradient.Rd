\name{gpsimMapWGradient}
\Rdversion{1.0}
\alias{gpsimMapWGradient}
\title{Compute the gradients of W with respect to the parameters of the k-th gene.}
\description{
  computes the gradients of W with respect to the parameters of the k-th gene given Gaussian process for use in a single input motif protein network.
}
\usage{
values <- gpsimMapWGradient(model, k)
}
\arguments{
  \item{model}{the model for which the log likelihood is computed.}
  \item{k}{the k-th gene.}
}
\value{
  \item{dWdB}{the gradients of W w.r.t the model paramter Bk.}
  \item{dWdD}{the gradients of W w.r.t the model paramter Dk.}
  \item{dWdS}{the gradients of W w.r.t the model paramter Sk.}
  \item{dWdalpha}{the gradients of W w.r.t the model paramter alpha_k.}
  \item{dWdgParam}{the gradients of W w.r.t the model paramter gamma_k}
  \item{dWdn}{the gradients of W w.r.t the general noise variance of the k-gene.}
}
\seealso{
\code{\link{gpsimMapCreate}, \link{gpsimMapLogLikelihood}, \link{gpsimMapGradient}, \link{gpsimMapFunctionalLogLikeGradients}}.
}
\examples{
## missing
}
\keyword{model}
