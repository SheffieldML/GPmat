\name{SCGoptim}
\Rdversion{1.0}
\alias{SCGoptim}
\alias{CGoptim}
\alias{optimiDefaultOptions}
\alias{modelOptimise}
\alias{cgpsimOptimise}

\title{Optimise the given function using (scaled) conjugate gradients.}
\description{
  Optimise the given function using (scaled) conjugate gradients.
}
\usage{
  options <- optimiDefaultOptions()
  newParams <- SCGoptim(x, fn, grad, options, ...)
  newParams <- CGoptim(x, fn, grad, options, ...)
  model <- modelOptimise(model, options, ...)
}
\arguments{
  \item{model}{the model to be optimised.}
  \item{x}{initial parameter values.}
  \item{fn}{objective function to minimise}
  \item{grad}{gradient function of the objective}
  \item{options}{options structure like one returned by
    \code{optimiDefaultOptions}. The fields are interpreted as\\
    option[1] : number of iterations\\
    option[2] : interval for the line search\\
    option[3] : tolerence for x to terminate the loop\\
    option[4] : tolerence for fn to terminate the loop\\
    option$display : option of showing the details of optimisaton
  }
  \item{...}{extra arguments to pass to fn and grad}
}
\value{
  \item{options}{an options structure}
  \item{newParams}{optimised parameter values}
  \item{model}{the optimised model.}
}
\seealso{
\code{\link{modelObjective}, \link{modelGradient}}
}
\examples{
## Not run to speed up package checks
# model <- GPLearn(..., dontOptimise=TRUE)
# options <- optimiDefaultOptions()
# model <- modelOptimise(model, options)
}
\keyword{model}
