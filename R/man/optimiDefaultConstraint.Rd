\name{optimiDefaultConstraint}
\Rdversion{1.0}
\alias{optimiDefaultConstraint}
\title{Returns function for parameter constraint.}
\description{
  returns the current default function for constraining a
  parameter.
}
\usage{
  optimiDefaultConstraint(constraint)
}
\arguments{
  \item{constraint}{the type of constraint you want to place on the parameter, options include 'positive' (gives an 'exp' constraint) and 'zeroone' (gives a 'sigmoid' constraint).}
}
\value{
  \item{val}{a list with two components: 'func' for the name of function
    used to apply the constraint, and 'hasArgs' for a boolean flag if
    the function requires additional arguments.}
}
\seealso{
\code{\link{expTransform}, \link{sigmoidTransform}}.
}
\examples{
  optimiDefaultConstraint('positive')
  optimiDefaultConstraint('bounded')
}
\keyword{model}
