\name{gpPosteriorSample}
\Rdversion{1.0}
\alias{gpPosteriorSample}
\title{Plot Samples from a GP Posterior.}
\description{
  Gaussian processes are non-parametric models. They are specified
  by their covariance function and a mean function. When combined
  with data observations a posterior Gaussian process is induced.
  This function samples from that posterior.
}
\usage{
  gpPosteriorSample(kernType, numSamps=10, params=NULL, lims=c(-3,3), path=getwd(), png=FALSE)
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
  \code{
    \link{gpOptions}, \link{kernCreate}, \link{kernCompute},
    \link{gaussSamp}, \link{zeroAxes}.
  }
}

