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
  
  \item{paramsList}{indices of parameteres to group together. The
  indices are provided in a list. Each element in the list contains
  a vector of indices of parameters that should be considered as one
  parameter. Each group of parameters in each cell should obviously be
  mutually exclusive.
  
  Alternatively, the specification may consist of strings, which are
  interpreted as regular expressions that are matched against the
  parameter names returned by \code{modelExtractParam} or
  \code{kernExtractParam}, as appropriate fot the current object.}
}
\value{
  \item{model}{the model with the parameters grouped together.}
}
\seealso{
\code{\link{modelExtractParam}, \link{modelExpandParam}, \link{modelGradient}}.
}
\examples{
# Create a multi kernel with two rbf blocks with bounded inverse widths
invWidthBounds <- c(0.5, 2)
kernType <- list(type="multi", comp=list())
for (i in 1:2)
  kernType$comp[[i]] <- list(type="parametric", realType="rbf",
                              options=list(isNormalised=TRUE,
                                inverseWidthBounds=invWidthBounds))
kern <- kernCreate(1, kernType)

# Tie the inverse with parameters of the component RBF kernels
kern <- modelTieParam(kern, list(tieWidth="inverseWidth"))
kernDisplay(kern)
}
\keyword{model}
