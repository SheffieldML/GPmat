\name{gpsimLoadBarencoPUMAData}
\Rdversion{1.0}
\alias{gpsimLoadBarencoPUMAData}
\title{Load in Martino Barenco's data as re-processed by mmgMOS.}
\description{
  loads in from the Excel spread sheets the data from the Barenco et al paper as processed by mmgMOS. OPTION defines the way of transforming the data from logged to un-logged space.
}
\usage{
values <- gpsimLoadBarencoPUMAData()
}
\arguments{
}
\value{
  \item{y}{the normalised expression levels.}
  \item{yvar}{the variance of the normalised expression levels.}
  \item{gene}{the gene names and Affymetrix array tags.}
  \item{times}{the times of the expression measurements.}
  \item{scale}{the scaling factor applied to normalise.}
  \item{rawExp}{the raw gene expresion level.}
  \item{rawVar}{the raw variance of the gene expression.}
}
\seealso{
\code{\link{demBarenco1}, \link{demBarencoMap1}}.
}
\examples{
## missing
}
\keyword{model}
