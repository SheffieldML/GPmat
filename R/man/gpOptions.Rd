\name{gpOptions}
\Rdversion{1.0}
\alias{gpOptions}
\alias{mlpOptions}
\title{Return default options for GP model.}
\description{
  returns the default options in a structure for a GP model.
}
\usage{
gpOptions(approx)
}
\arguments{
  \item{approx}{approximation type, either 'ftc' (no approximation), 'dtcvar' (variational sparse approximation) 'dtc' (deterministic training conditional), 'fitc' (fully independent training conditional) or 'pitc' (partially independent training conditional).}
}
\value{
  \item{options}{structure containing the default options for the given approximation type.}
}
\seealso{
\code{\link{gpCreate}}.
}
\examples{
## missing
}
\keyword{model}
