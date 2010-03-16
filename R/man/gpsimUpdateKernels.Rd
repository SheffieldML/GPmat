\name{gpsimUpdateKernels}
\Rdversion{1.0}
\alias{gpsimUpdateKernels}
\title{Updates the kernel representations in the GPSIM structure.}
\description{
  updates any representations of the kernel in the model structure, such as invK, logDetK or K.
}
\usage{
model <- gpsimUpdateKernels(model)
}
\arguments{
  \item{model}{the model structure for which kernels are being updated.}
}
\value{
  \item{model}{the model structure with the kernels updated.}
}
\seealso{
\code{\link{modelUpdateProcesses}}.
}
\examples{
## missing
}
\keyword{model}
