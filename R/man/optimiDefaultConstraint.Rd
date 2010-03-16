\name{optimiDefaultConstraint}
\Rdversion{1.0}
\alias{optimiDefaultConstraint}
\title{Returns function for parameter constraint.}
\description{
  returns the current default function for constraining a
  parameter.
}
\usage{
str <- optimiDefaultConstraint(constraint)
}
\arguments{
  \item{constraint}{the type of constraint you want to place on the parameter, options include 'positive' (gives an 'exp' constraint) and 'zeroone' (gives a 'sigmoid' constraint).}
}
\value{
  \item{str}{the type of function used to apply the constraint from the 'optimi' toolbox.}
}
\seealso{
\code{\link{expTransform}, \link{sigmoidTransform}}.
}
\examples{
## missing
}
\keyword{model}
