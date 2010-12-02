\name{modelGradientCheck}
\Rdversion{1.0}
\alias{modelGradientCheck}
\title{Check gradients of given model.}
\description{
  checks the supplied gradient function and the supplied objective function to ensure that the numerical gradients (as computed with the objective function) match the analytically computed gradients.
}
\usage{
modelGradientCheck(model, P1, P2, P3 ...)
}
\arguments{
  \item{model}{the model for which gradients are to be checked.}
  \item{P1, P2, P3 ...}{additional arguments that are passed to the objective and gradient functions (after the parameter vector which is always assumed to be the first argument passed).}
}
\value{
}
\seealso{
\code{\link{modelObjective}, \link{modelGradient}, \link{modelCreate}}.
}
\examples{
## missing
}
\keyword{model}
