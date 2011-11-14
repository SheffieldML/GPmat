\name{modelDisplay}
\Rdversion{1.0}
\alias{modelDisplay}
\alias{gpdisimDisplay}
\alias{gpsimDisplay}
\alias{cmpndKernDisplay}
\alias{disimKernDisplay}
\alias{kernDisplay}
%\alias{mlpKernDisplay}
\alias{multiKernDisplay}
\alias{rbfKernDisplay}
\alias{simKernDisplay}
%\alias{translateKernDisplay}
\alias{whiteKernDisplay}
\title{Display a model.}
\description{
  displays the parameters of the model/kernel and the model/kernel type to
  the console.
}
\usage{
  modelDisplay(model, ...)
}
\arguments{
  \item{model}{the model/kernel structure to be displayed.}
  \item{...}{(optional) indent level for the display.}
}
\seealso{
\code{\link{modelExtractParam}}
}
\examples{
  # Load a mmgmos preprocessed fragment of the Drosophila developmental
  # time series
  data(drosophila_gpsim_fragment)

  # The probe identifier for TF 'twi'
  twi <- "143396_at"
  # The probe identifier for the target gene
  targetProbe <- "152715_at"

  # Create the model, but do not optimise
  model <- GPLearn(drosophila_gpsim_fragment,
                   TF=twi, targets=targetProbe,
                   useGpdisim=TRUE, quiet=TRUE,
                   dontOptimise=TRUE)

  # Display the initial model
  modelDisplay(model)
}
\keyword{model}
