\name{gpCreate}
\Rdversion{1.0}
\alias{gpCreate}
\title{Create a GP model with inducing varibles/pseudo-inputs.}
\description{
  Creates a Gaussian process model structure with default parameter settings as specified by the options vector.
}
\usage{
model <- gpCreate(q, d, X, y, options)
}
\arguments{
  \item{q}{input data dimension.}
  \item{d}{the number of processes (i.e. output data dimension).}
  \item{X}{the input data matrix.}
  \item{y}{the target (output) data.}
  \item{options}{options structure as defined by gpOptions.R.}
}
\value{
  \item{model}{model structure containing the Gaussian process.}
}
\seealso{
\code{\link{gpOptions}.
}
}
\examples{
## missing
}
\keyword{model}
