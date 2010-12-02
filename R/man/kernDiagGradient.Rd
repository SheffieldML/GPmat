\name{kernDiagGradient}
\Rdversion{1.0}
\alias{kernDiagGradient}
\title{Compute the gradient of the kernel's parameters for the diagonal.}
\description{
  computes the gradient of functions of the diagonal of the kernel matrix with respect to the parameters of the kernel. The parameters' gradients are returned in the order given by the kernExtractParam command.
}
\usage{
g <- kernDiagGradient(kern, x, factors)
}
\arguments{
  \item{kern}{the kernel structure for which the gradients are computed.}
  \item{x}{the input data for which the gradient is being computed.}
  \item{factors}{partial derivatives of the function of interest with respect to the diagonal elements of the kernel matrix.}
}
\value{
  \item{g}{gradients of the relevant function with respect to each of the parameters. Ordering should match the ordering given in kernExtractParam.}
}
\seealso{
\code{\link{kernDiagGradient}, \link{kernExtractParam}, \link{kernGradient}}.
}
\examples{
## missing
}
\keyword{model}
