\name{expTransform}
\Rdversion{1.0}
\alias{expTransform}
\title{Constrains a parameter to be positive through exponentiation.}
\description{
  contains commands to constrain parameters to be positive via exponentiation.
}
\usage{
expTransform(x, y, transform)
}
\arguments{
  \item{x}{input argument.}
  \item{y}{return argument.}
  \item{transform}{type of transform, 'atox' maps a value into the transformed space (i.e. makes it positive). 'xtoa' maps the parameter back from transformed space to the original space. 'gradfact' gives the factor needed to correct gradients with respect to the transformed parameter. }
}
\value{
}
\seealso{
\code{\link{negLogLogitTransform}, \link{sigmoidTransform}}.
}
\examples{
## missing
}
\keyword{model}
