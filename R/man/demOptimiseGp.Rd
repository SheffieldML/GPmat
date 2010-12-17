\name{demOptimiseGp}
\Rdversion{1.0}
\alias{demOptimiseGp}
\title{Gaussian Process Optimisation Demo}
\description{
  Shows a series of plots of a Gaussian process with different length
  scales fitted to six data points. For each plot there is a corresponding
  plot of the log likelihood. The log likelihood peaks for a length scale
  close to 1. This was the length scale used to generate the data.
}
\usage{
  demOptimiseGp(path=getwd(), filename='demOptimiseGp', png=FALSE, gif=FALSE)
}
\arguments{
  \item{path}{path where the plot images are saved.}
  \item{filename}{name of saved images.}
  \item{png}{save image as png.}
  \item{gif}{save series of images as animated gif.}
}
\seealso{
  \code{
    \link{gpOptions}, \link{kernCreate}, \link{kernCompute},
    \link{gaussSamp}, \link{kernDiagCompute}, \link{gpCreate},
    \link{gpExpandParam}, \link{gpLogLikelihood}, \link{gpPosteriorMeanVar}, \link{gpPlot}.
  }
}
\keyword{model}

