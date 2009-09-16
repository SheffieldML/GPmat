\name{gpsimMapUpdateKernels}
\Rdversion{1.0}
\alias{gpsimMapUpdateKernels}
\title{Updates the kernel representations in the GPSIMMAP structure.}
\description{
  updates any representations of the kernel in the model structure, such as invK, logDetK or K.
}
\usage{
values <- gpsimMapUpdateKernels(model)
}
\arguments{
  \item{model}{the model structure for which kernels are being updated.}
}
\value{
  \item{model}{the model structure with the kernels updated.}
  \item{ll}{the log likelihood of the model after update.}
}
\seealso{
\code{\link{}}.
}
\examples{
## missing
}
\keyword{model}
