\name{gpOptimise}
\Rdversion{1.0}
\alias{gpOptimise}
\title{Optimise the inducing variable based kernel.}
\description{
  optimises the Gaussian process  model for a given number of iterations.
}
\usage{
model <- gpOptimise(model, display, iters)
}
\arguments{
  \item{model}{the model to be optimised.}
  \item{display}{whether or not to display while optimisation proceeds, set to 2 for the most verbose and 0 for the least verbose.}
  \item{iters}{number of iterations for the optimisation.}
}
\value{
  \item{model}{the optimised model.}
}
\seealso{
\code{
  \link{gpCreate}, \link{gpGradient}, \link{gpObjective}}.
}
\examples{
## missing
}
\keyword{model}
