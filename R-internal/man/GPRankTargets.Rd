\name{GPRankTargets}
\Rdversion{1.0}
\alias{GPRankTargets}
\title{Ranking possible target genes by forming optimized models}
\description{
  Ranks possible target genes by forming optimized models with a fixed
  transcription factor, a set of known target genes and targets to be
  tested. The transcription factor and the known targets are always
  included in the models while the tested targets are tested by
  including them in the models one at a time. The function determines
  itself whether to use GPsim or GPdisim based on the input
  arguments. Finally, the formed models are sorted based on the
  logarithmic value of their marginal likelihood.
}
\usage{
scores <- GPRankTargets(preprocData, TF = NULL, knownTargets = NULL,
testTargets = NULL, filterLimit = 1.8, returnModels = FALSE, options =
NULL, scoreSaveFile = NULL)
}
\arguments{
  \item{preprocData}{The preprocessed data to be used.}
  \item{TF}{The transcription factor present in all models.}
  \item{knownTargets}{The target genes present in all models.}
  \item{testTargets}{Target genes that are tested by including them in
    the models one at a time. Can be names of genes, or a set of indices
    in preprocData.}
  \item{filterLimit}{Genes with an average expression z-score above
    this figure are accepted after filtering. If this value is 0, all
    genes will be accepted.}
  \item{returnModels}{A logical value determining whether the function
    returns the calculated models.}
  \item{options}{A list of additional arguments to pass to GPLearn.}
  \item{scoreSaveFile}{Name of file to save the scores to after processing
    each gene.}
}

\details{The models are formed by calling the function formModel, which
  in turn calls GPLearn. If there is no value given to the transcription
  factor, GPSIM is used. Otherwise, GPDISIM is used. GPSIM needs some
  known targets. If known
  targets are given, a model is first created with only the transcription
  factor and the known targets. The parameters extracted from this model
  are used as initial parameters of the models with test targets.
}

\value{
  The function returns a scoreList
  containing the genes, parameters and log-likelihoods of the models
  If returnModels is true,
  the function returns a list of the calculated models.
}
\author{
  Jonatan Ropponen, Antti Honkela
}
\seealso{
  \code{\link{GPLearn}, \link{GPRankTFs}, \linkS4class{scoreList}, \link{generateModels}}.
}
\examples{
  ## Load a mmgmos preprocessed fragment of the Drosophila developmental
  ## time series
  data(drosophila_gpsim_fragment)

  ## Get the target probe names
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
}
\keyword{model}
