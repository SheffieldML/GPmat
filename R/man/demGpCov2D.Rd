\name{demGpCov2D}
\Rdversion{1.0}
\alias{demGpCov2D}
\title{Gaussian Process 2D Covariance Demo}
\description{
  Gives the joint distribution for f1 and f2, two values of a function sampled from a Gaussian process.
  The plots show the joint distributions as well as the conditional for f2 given f1.
}
\usage{
  demGpCov2D(ind=c(1,2), path = getwd(), filename = paste('demGpCov2D', ind[1],'_', ind[2], sep=''), png=FALSE, gif=FALSE)
}
\arguments{
  \item{ind}{indices of function values.}
  \item{path}{path where the plot images are saved.}
  \item{filename}{name of saved images.}
  \item{png}{save image as png.}
  \item{gif}{save series of images as animated gif.}
}
\seealso{
  \code{\link{demGpSample}, \link{basePlot}, \link{zeroAxes}.}
}
\keyword{model}

