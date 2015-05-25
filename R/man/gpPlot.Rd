\name{gpPlot}
\Rdversion{1.0}
\alias{gpPlot}
\title{Gaussian Process Plotter}
\description{
  Plots the GP mean and variance.
}
\usage{
  gpPlot(model, Xstar, mu, S, simpose=NULL,xlim=NULL,ylim=NULL,xlab='',ylab='',col='blue',title='')
}
\arguments{
  \item{model}{the model structure for which GP mean and variance are to be plotted.}
  \item{Xstar}{the input positions for which the mean and variance will be plotted.}
  \item{mu}{the precomputed GP posterior mean vector.}
  \item{S}{the precomputed GP posterior variance vector.}
  \item{simpose}{vector of datapoints to be superimposed on the plot with added white noise.}
  \item{xlim}{x-axis plotting limits.}
  \item{ylim}{y-axis plotting limits.}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{col}{color for plotting the GP mean and variance.}
  \item{title}{plot title.}
}
\seealso{
  \code{
    \link{gpPosteriorMeanVar}, \link{polygon}, \link{zeroAxes}.
  }
}
\keyword{model}

