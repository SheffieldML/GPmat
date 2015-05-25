\name{GPPlot}
\Rdversion{1.0}
\alias{GPPlot}
\title{Plot GP(DI)SIM models}
\description{
  Plots GP(DI)SIM models.
}
\usage{
GPPlot(data, savepath = '', nameMapping = NULL, predt = NULL,
       fileOutput=FALSE, plotTime=NULL)
}
\arguments{
  \item{data}{The model to plot as returned by GPLearn.}
  \item{savepath}{The location in the file system where the images are saved.}
  \item{nameMapping}{The annotation used for mapping the names of the
    genes for the figures.}
  \item{predt}{The set of time points to use in plotting
    (default: the time interval covering the data).}
  \item{fileOutput}{Is the plot being saved to a file? If yes, do not
    open new interactive devices for each plot.}
  \item{plotTime}{The times of observations to use in the plot.
    Should usually not be changed!}
}
\details{
  The function plots
  the fitted expression level of the transcription factor (if applicable),
  the inferred activity of the transcription factor,
  and the fitted expression level of the target(s).
}
\author{Antti Honkela}
\seealso{
  \code{\link{GPLearn}}.
}
\examples{\dontrun{
  # Load a mmgmos preprocessed fragment of the Drosophila developmental
  # time series
  data(drosophila_gpsim_fragment)

  # Get the target probe names
  library(annotate)
  aliasMapping <- getAnnMap("ALIAS2PROBE",
                    annotation(drosophila_gpsim_fragment))
  twi <- get('twi', env=aliasMapping)
  fbgnMapping <- getAnnMap("FLYBASE2PROBE",
                   annotation(drosophila_gpsim_fragment))
  targetProbe <- get('FBgn0035257', env=fbgnMapping)

  # Learn the model
  model <- GPLearn(drosophila_gpsim_fragment,
                   TF=twi, targets=targetProbe,
                   useGpdisim=TRUE, quiet=TRUE)

  # Plot it
  GPPlot(model, nameMapping=getAnnMap("FLYBASE",
         annotation(drosophila_gpsim_fragment)))
}}
\keyword{model}
