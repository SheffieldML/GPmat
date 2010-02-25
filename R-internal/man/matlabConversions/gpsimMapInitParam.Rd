\name{gpsimMapInitParam}
\Rdversion{1.0}
\alias{gpsimMapInitParam}
\title{Creates a set of options for a GPSIMMAP model as the initial parameters.}
\description{
  returns a default initialisation for a GPSIMMAP model.
}
\usage{
options <- gpsimMapInitParam(y, yvar, options)
}
\arguments{
  \item{y}{the values of each gene at the different time points.}
  \item{yvar}{the variances of each gene at the different time points.}
  \item{options}{the values of each gene at the different time points.}
}
\value{
  \item{options}{the options structure. }
}
\seealso{
\code{\link{gpsimMapCreate}}.
}
\examples{
## missing
}
\keyword{model}
