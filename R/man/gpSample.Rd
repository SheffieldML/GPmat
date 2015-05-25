\name{gpSample}
\Rdversion{1.0}
\alias{gpSample}
\title{Plot Samples from a GP.}
\description{
  creates a plot of samples from a kernel with the given parameters and variance.
}
\usage{
gpSample(kernType, numSamps=10, params=NULL, lims=c(-3,3), path=getwd(), png=FALSE)
}
\arguments{
  \item{kernType}{the type of kernel to sample from.}
  \item{numSamps}{the number of samples to take.}
  \item{params}{parameter vector for the kernel.}
  \item{lims}{limits of the x axis.}
  \item{path}{path where the plot images are saved.}
  \item{png}{save image as png.}
}
\seealso{
  \code{\link{gpOptions}.
}
}

