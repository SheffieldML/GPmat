\name{gpExpandParam}
\Rdversion{1.0}
\alias{gpExpandParam}
\title{Expand a parameter vector into a GP model.}
\description{
  takes the given vector of parameters and places them in the model structure, it then updates any stored representations that are dependent on those parameters, for example kernel matrices etc..
}
\usage{
gpExpandParam(model, params)
}
\arguments{
  \item{model}{the model structure for which parameters are to be updated.}
  \item{params}{a vector of parameters for placing in the model structure.}
}
\value{
  \item{model}{a returned model structure containing the updated parameters. }
}
\seealso{
\code{\link{gpCreate}, \link{gpExtractParam}, \link{modelExtractParam}, \link{gpUpdateKernels}}.
}
\examples{
## missing
}
\keyword{model}
