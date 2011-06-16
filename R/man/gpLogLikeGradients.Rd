\name{gpLogLikeGradients}
\Rdversion{1.0}
\alias{gpLogLikeGradients}
\alias{localCovarianceGradients}
\alias{localSCovarianceGradients}
\title{Compute the gradients for the parameters and X.}
\description{
  computes the gradients of the Gaussian process log likelihood with respect
  to the model parameters (and optionally, as above with respect to inducing
  variables and input data) given the target data, input data and inducing
  variable locations.
}
\usage{
gpLogLikeGradients(model, X=model$X, M, X_u, gX_u.return=FALSE, gX.return=FALSE, g_beta.return=FALSE)
}
\arguments{
  \item{model}{the model structure for which gradients are computed.}
  \item{X}{the input data locations for which gradients are computed.}
  \item{M}{the scaled and bias removed target data for which the gradients are computed.}
  \item{X_u}{the inducing variable locations for which gradients are computed.}
  \item{gX_u.return}{(logical) return the gradient of the log likelihood with respect to the inducing variables. If inducing variables aren't being used this returns zero.}
  \item{gX.return}{(logical) return the gradient of the log likelihood with respect to the input data locations.}
  \item{g_beta.return}{(logical) to return the gradient of the log likelihood with respect to beta.}
}
\value{
  \item{gParam}{contains the gradient of the log likelihood with respect to the model parameters (including any gradients with respect to beta).}
}
\seealso{
\code{
  \link{gpLogLikelihood}}.
}
\examples{
## missing
}
\keyword{model}
