\name{demRegression}
\Rdversion{1.0}
\alias{demRegression}
\title{Gaussian Process Regression Demo}
\description{
  The regression demo very much follows the format of the interpolation demo.
  Here the difference is that the data is sampled with noise. Fitting a model
  with noise means that the regression will not necessarily pass right through
  each data point.
}
\usage{
  demRegression(path=getwd(), filename='demRegression')
}
\arguments{
  \item{path}{path where the plot images are saved.}
  \item{filename}{name of saved images.}
}
\seealso{
  \code{
    \link{gpOptions}, \link{kernCreate}, \link{kernCompute},
    \link{gaussSamp}, \link{kernDiagCompute}, \link{gpCreate},
    \link{gpPlot}, \link{demInterpolation}.
  }
}
\keyword{model}

