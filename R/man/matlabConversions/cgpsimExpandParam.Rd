\name{cgpsimExpandParam}
\Rdversion{1.0}
\alias{cgpsimExpandParam}
\title{Expand params into model structure.}
\description{
  this is just a wrapper for learning three GP models simultaneously.
}
\usage{
model <- cgpsimExpandParam(model, params)
}
\arguments{
  \item{model}{input model from which parameters will be extracted.}
  \item{params}{parameters to place in the model.}
}
\value{
  \item{model}{model filled with the given parameters.}
}
\seealso{
\code{\link{gpsimExpandParam}}.
}
\examples{
## missing
}
\keyword{model}
