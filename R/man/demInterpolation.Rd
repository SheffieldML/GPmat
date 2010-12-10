\name{demInterpolation}
\Rdversion{1.0}
\alias{demInterpolation}
\title{Gaussian Process Interpolation Demo}
\description{
  Plots, consecutively, an increasing number of data points, followed
  by an interpolated fit through the data points using a Gaussian process.
  This is a noiseless system, and the data is sampled from a GP with a
  known covariance function. The curve is then recovered with minimal
  uncertainty after only nine data points are included.
}
\usage{
  demInterpolation(path=getwd(), filename='demInterpolation')
}
\arguments{
  \item{path}{path where the plot images are saved.}
  \item{filename}{name of saved images.}
}
\seealso{
  \code{
    \link{gpOptions}, \link{kernCreate}, \link{kernCompute},
    \link{gaussSamp}, \link{kernDiagCompute}, \link{gpCreate}, \link{gpPlot}.
  }
}
\keyword{model}

