\name{optimiDefaultConstraint}
\Rdversion{1.0}
\alias{optimiDefaultConstraint}
\title{Returns function for parameter constraint.}
\description{
  returns the current default function for constraining a parameter. Formerly (up to version KERN 0.163) this was 'negLogLogit' for positive constraints, as this keeps things roughly linear in the positive half space, however, it is more standard to use 'exp' (i.e. optimise in the log space). This function allows you to change the option globally. The defaults are given below.
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
\code{\link{expTransform}, \link{sigmoidTransform}, \link{linearTransform}, \link{negLogLogitTransform}}.
}
\examples{
## missing
}
\keyword{model}
