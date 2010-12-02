\name{gpUpdateKernels}
\Rdversion{1.0}
\alias{gpUpdateKernels}
\title{Update the kernels that are needed.}
\description{
  updates any representations of the kernel in the model structure, such as invK, logDetK or K.
}
\usage{
values <- gpUpdateKernels(model, model, X, X_u)
}
\arguments{
  \item{model}{the model structure for which kernels are being updated.}
  \item{model}{the model structure for which kernels are being updated.}
  \item{X}{the input locations for update of kernels.}
  \item{X_u}{the inducing input locations.}
}
\value{
  \item{model}{the model structure with the kernels updated.}
  \item{model}{the model structure with the kernels updated.}
}
\seealso{
\code{\link{gpExpandParam}, \link{gpCreate}}.
}
\examples{
## missing
}
\keyword{model}
