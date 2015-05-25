\name{demAutoOptimiseGp}
\Rdversion{1.0}
\alias{demAutoOptimiseGp}
\title{Gaussian Process Optimisation Demo}
\description{
  Shows that by varying the length scale, an artificial data set has
  different likelihoods, yet there is an optimum for which the likelihood is
  maximised. This demo is similar to demOptimiseGp, only here, it is
  demonstrated how the length scale hyperparameter is optimised automatically
  through SCG (scaled conjugate gradients) numerical optimisation.
  Run multiple times to understand the effect of optimisation on randomly generated datasets.
}
\usage{
  demAutoOptimiseGp(path=getwd(), filename='demAutoOptimiseGp', png=FALSE, gif=FALSE)
}
\arguments{
  \item{path}{path where the plot images are saved.}
  \item{filename}{name of saved images.}
  \item{png}{save image as png.}
  \item{gif}{save series of images as animated gif.}
}
\seealso{
  \code{
    \link{gpOptions}, \link{kernCreate},
    \link{gaussSamp}, \link{gpCreate},
    \link{gpExpandParam}, \link{gpLogLikelihood}, \link{gpPosteriorMeanVar},
    \link{gpOptimise}, \link{gpPlot}.
  }
}
\keyword{model}

