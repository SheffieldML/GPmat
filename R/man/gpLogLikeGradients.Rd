\name{gpLogLikeGradients}
\Rdversion{1.0}
\alias{gpLogLikeGradients}
\title{Compute the gradients for the parameters and X.}
\description{
  computes the gradients of the Gaussian process log likelihood with respect to the model parameters (and optionally, as above with respect to inducing variables and input data) given the target data, input data and inducing variable locations.
}
\usage{
values <- gpLogLikeGradients(, , , , , ', ', )
}
\arguments{
  \item{}{model : the model structure for which gradients are computed.}
  \item{'}{X : the input data locations for which gradients are computed.}
  \item{'}{M : the scaled and bias removed target data for which the gradients are computed.}
  \item{'}{X_U : the inducing variable locations for which gradients are computed.}
}
\value{
  \item{gParam}{the gradient of the log likelihood with respect to the model parameters.}
  \item{gX_u}{the gradient of the log likelihood with respect to the inducing variables. If inducing variables aren't being used this returns zero.}
  \item{gParam}{the gradient of the log likelihood with respect to the model parameters (including any gradients with respect to beta).}
  \item{gX_u}{the gradient of the log likelihood with respect to the inducing variables. If inducing variables aren't being used this returns zero.}
  \item{gX}{the gradient of the log likelihood with respect to the input data locations.}
  \item{gbeta}{the gradient of the log likelihood with respect to beta.}
  \item{gParam}{the gradient of the log likelihood with respect to the model parameters.}
}
\seealso{
\code{\link{gpLogLikelihood}, \link{modelLogLikeGradients}, \link{fgplvmLogLikeGradients}}.
}
\examples{
## missing
}
\keyword{model}
