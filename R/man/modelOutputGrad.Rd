\name{modelOutputGrad}
\Rdversion{1.0}
\alias{modelOutputGrad}
\title{Compute derivatives with respect to params of model outputs.}
\description{
  Compute derivatives with respect to params of model outputs.
}
\usage{
  modelOutputGrad(model, X, dim)
}
\arguments{
  \item{model}{the model structure for which gradients are computed.}
  \item{X}{input locations where gradients are to be computed.}
  \item{dim}{the dimension of the model for which gradients are required.}
}
\details{
  \code{g <- modelOutputGrad(model, X)}
  gives the gradients of the outputs from the model with respect to the parameters for a given set of inputs.

  \code{g <- modelOutputGrad(model, X, dim)}
  gives the gradients of the outputs from the model with respect to the parameters for a given set of inputs.

}
\value{
  \item{g}{gradients of the model output with respect to the model parameters for the given input locations. The size of the returned matrix is of dimension number of data x number of parameters x number of model outputs (which maintains compatability with NETLAB).}
}
\seealso{
\code{
  \link{modelLogLikelihood}}.
}
\examples{
## missing
}
\keyword{model}
