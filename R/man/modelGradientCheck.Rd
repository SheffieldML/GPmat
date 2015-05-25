\name{modelGradientCheck}
\Rdversion{1.0}
\alias{modelGradientCheck}
\title{Check gradients of given model.}
\description{
  checks the supplied gradient function and the supplied objective function to ensure that the numerical gradients (as computed with the objective function) match the analytically computed gradients.
}
\usage{
  modelGradientCheck(model, ...)
}
\arguments{
  \item{model}{the model for which gradients are to be checked.}
  \item{...}{additional arguments that are passed to the objective and gradient functions (after the parameter vector which is always assumed to be the first argument passed).}
}
\seealso{
\code{
  \link{modelObjective}, \link{modelGradient}}.
}
\keyword{model}
