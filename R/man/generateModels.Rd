\name{generateModels}
\Rdversion{1.0}
\alias{generateModels}
\title{Generating models with the given data}
\description{
  'generateModels' recreates models based on the parameters
  stored in a scoreList.
}
\usage{
models <- generateModels(preprocData, scores)
}
\arguments{
  \item{preprocData}{The preprocessed data to be used.}
  \item{scores}{A scoreList object containing data of the models to be
    generated.}
}
%\details{
%  'generateModels' first creates a model with 'gpsimCreate' or 'gpdisimCreate'. Then 'gpsimExpandParam' or 'gpdisimExpandParam' is used to generate the original models without having to optimize them again.
%}
\value{
  'generateModels' returns a list of the generated models.
}
\author{Antti Honkela, Jonatan Ropponen}
\seealso{
  \code{\link{GPLearn}, \link{GPRankTargets}, \link{GPRankTFs}, \linkS4class{scoreList}}.
}
\examples{\dontrun{
  # Load a mmgmos preprocessed fragment of the Drosophila developmental
  # time series
  data(drosophila_gpsim_fragment)

  # Get the target probe names
  targets <- c('FBgn0003486', 'FBgn0033188', 'FBgn0035257')
  library(annotate)
  aliasMapping <- getAnnMap("ALIAS2PROBE",
                    annotation(drosophila_gpsim_fragment))
  twi <- get('twi', env=aliasMapping)
  fbgnMapping <- getAnnMap("FLYBASE2PROBE",
                   annotation(drosophila_gpsim_fragment))
  targetProbes <- mget(targets, env=fbgnMapping)

  scores <- GPRankTargets(drosophila_gpsim_fragment, TF=twi,
                          testTargets=targetProbes,
                          options=list(quiet=TRUE),
                          filterLimit=1.8)

  models <- generateModels(drosophila_gpsim_fragment, scores)
}}
\keyword{model}
