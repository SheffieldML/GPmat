\name{noiseCreate}
\Rdversion{1.0}
\alias{noiseCreate}
\title{Initialise a noise structure.}
\description{
  takes a noise type and a target vector and initialises a noise structure from it. The parameters of the noise structure are the set by calling noiseParamInit.
}
\usage{
noiseCreate(noiseType, y)
}
\arguments{
  \item{noiseType}{the type of noise to be created (e.g. 'gaussian', 'probit', 'ncnm').}
  \item{y}{the target vector.}
}
\value{
}
\seealso{
\code{\link{noiseParamInit}}.
}
\examples{
## missing
}
\keyword{model}
