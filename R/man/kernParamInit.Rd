\name{kernParamInit}
\Rdversion{1.0}
\alias{kernParamInit}
\alias{cmpndKernParamInit}
\alias{disimKernParamInit}
\alias{mlpKernParamInit}
\alias{multiKernParamInit}
\alias{rbfKernParamInit}
\alias{simKernParamInit}
\alias{translateKernParamInit}
\alias{whiteKernParamInit}

\title{Kernel parameter initialisation.}
\description{
  initialises the parameters of a kernel.
}
\usage{
kern <- kernParamInit(kern)
}
\arguments{
  \item{kern}{the kernel structure for which the parameters will be initialised.}
}
\value{
  \item{kern}{the kernel structure with the parameters initialised.}
}
\seealso{
\code{\link{kernCreate}}.
}
\examples{
## missing
}
\keyword{model}
