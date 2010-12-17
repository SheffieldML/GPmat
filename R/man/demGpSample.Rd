\name{demGpSample}
\Rdversion{1.0}
\alias{demGpSample}
\title{Gaussian Process Sampling Demo}
\description{
  This example shows how points which look like they come from a function
  to be sampled from a Gaussian distribution. The sample is 25 dimensional
  and is from a Gaussian with a particular covariance.
}
\usage{
  demGpSample(bw=FALSE, path=getwd(), filename='gpSample', png=FALSE)
}
\arguments{
  \item{bw}{greyscale plots.}
  \item{path}{path where the plot images are saved.}
  \item{filename}{name of saved images.}
  \item{png}{save image as png.}
}
\seealso{
  \code{\link{kernCreate}, \link{kernCompute}, \link{gaussSamp}, \link{demGpSample}, \link{basePlot}, \link{zeroAxes}.}
}
\keyword{model}

