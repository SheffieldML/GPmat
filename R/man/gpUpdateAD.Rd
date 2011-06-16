\name{gpUpdateAD}
\Rdversion{1.0}
\alias{gpUpdateAD}
\title{Update the representations of A and D associated with the model.}
\description{
  updates the representations of A and D in the model when called by gpUpdateKernels.
}
\usage{
gpUpdateAD(model, X)
}
\arguments{
  \item{model}{the model for which the representations are being updated.}
  \item{X}{the X values for which the representations are being computed.}
}
\value{
  \item{model}{the model with the A and D representations updated.}
}
\seealso{
\code{\link{gpUpdateKernels}, \link{gpExpandParam}}.
}
\examples{
## missing
}
\keyword{model}
