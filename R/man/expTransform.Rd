\name{expTransform}
\Rdversion{1.0}
\alias{expTransform}
\alias{sigmoidTransform}
\alias{boundedTransform}
\title{Constrains a parameter.}
\description{
  contains commands to constrain parameters to be positive via
  exponentiation or within a fixed interval via the sigmoid function.
}
\usage{
expTransform(x, transform)
sigmoidTransform(x, transform)
boundedTransform(x, transform, bounds)
}
\arguments{
  \item{x}{input argument.}
  \item{transform}{type of transform, 'atox' maps a value into the
  transformed space (i.e. makes it positive). 'xtoa' maps the parameter
  back from transformed space to the original space. 'gradfact' gives
  the factor needed to correct gradients with respect to the transformed
  parameter.}
  \item{bounds}{a 2-vector of bounds of allowed values in boundedTransform}
}
\value{
  Return value as selected by \code{tranform}
}
\seealso{
\code{\link{modelOptimise}}
}
\examples{
  # Transform unconstrained parameter -4 to a positive value
  expTransform(-4, 'atox')

  # Transform a bounded parameter in (1,3) to an unconstrained one
  boundedTransform(2, 'xtoa', c(1, 3))
}
\keyword{model}
