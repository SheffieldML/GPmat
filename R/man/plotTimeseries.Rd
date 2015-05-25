\name{plotTimeseries}
\Rdversion{1.0}
\alias{plotTimeseries}
\title{Plot ExpressionTimeSeries data}
\description{
  Plots ExpressionTimeSeries data.
}
\usage{
plotTimeseries(data, nameMapping = NULL)
}
\arguments{
  \item{data}{An ExpressionTimeSeries object.}
  \item{nameMapping}{The annotation used for mapping the names of the
    genes for the figures.  By default, the SYMBOL annotation for the
    array is used, if available.}
}
\details{
  The function plots the expression levels from an ExpressionTimeSeries
  object and the associated standard deviations.  If the object
  includes multiple time series, they will be plotted in the same
  figure, but slightly shifted.
}
\author{Antti Honkela}
\seealso{
  \code{\link{processData}}.
}
\examples{
  # Load a mmgmos preprocessed fragment of the Drosophila developmental
  # time series
  data(drosophila_gpsim_fragment)

  # Plot the first two genes
  plotTimeseries(drosophila_gpsim_fragment[1:2,])
}
\keyword{model}
