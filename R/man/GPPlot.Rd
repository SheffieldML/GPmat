\name{GPPlot}
\Rdversion{1.0}
\alias{GPPlot}
\title{Plot GP(DI)SIM models}
\description{
  Plots GP(DI)SIM models.
}
\usage{
GPPlot(data, savepath = '', nameMapping = NULL, predt = NULL)
}
\arguments{
  \item{model}{The model to plot as returned by GPLearn.}
  \item{savepath}{The location in the file system where the images are saved.}
  \item{nameMapping}{The annotation used for mapping the names of the
    genes for the figures.}
  \item{predt}{The set of time points to use in plotting
    (default: the time interval covering the data).}
}
\details{
  The function plots
  the fitted expression level of the transcription factor (if applicable),
  the inferred activity of the transcription factor,
  and the fitted expression level of the target(s).
}
\author{
  Antti Honkela
}
\seealso{
  \code{\link{GPLearn}}.
}
\examples{
## missing
}
\keyword{model}
