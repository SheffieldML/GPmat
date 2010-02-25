\name{modelTieParam}
\Rdversion{1.0}
\alias{modelTieParam}
\title{Tie parameters of a model together.}
\description{
  groups of parameters of a model to be seen as one parameter during optimisation of the model.
}
\usage{
model <- modelTieParam(model, paramsList)
}
\arguments{
  \item{model}{the model for which parameters are being tied together.}
  \item{paramsList}{indices of parameteres to group together. The indices are provided in a cell array. Each cell in the array contains a vector of indices of parameters that should be considered as one parameter. Each group of parameters in each cell should obviously be mutually exclusive.}
}
\value{
  \item{model}{the model with the parameters grouped together.}
}
\seealso{
\code{\link{modelExtractParam}, \link{modeExpandParam}, \link{modelLogLikeGradients}}.
}
\examples{
## missing
}
\keyword{model}
