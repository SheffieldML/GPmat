\name{GPRankTargets}
\Rdversion{1.0}
\alias{GPRankTargets}
\alias{GPRankTFs}
\title{Ranking possible target genes or regulators}
\description{
  \code{GPRankTargets} ranks possible target genes by forming
  optimized models with a fixed
  transcription factor, a set of known target genes and targets to be
  tested. The transcription factor and the known targets are always
  included in the models while the tested targets are tested by
  including them in the models one at a time. The function determines
  itself whether to use GPSIM or GPDISIM based on the input
  arguments.
}
\usage{
scores <- GPRankTargets(preprocData, TF = NULL, knownTargets = NULL,
                        testTargets = NULL, filterLimit = 1.8,
                        returnModels = FALSE, options = NULL,
                        scoreSaveFile = NULL)
scores <- GPRankTFs(preprocData, TFs, targets,
                    filterLimit = 1.8, returnModels = FALSE, options = NULL,
                    scoreSaveFile = NULL)
}
\arguments{
  \item{preprocData}{The preprocessed data to be used.}
  \item{TF}{The transcription factor (TF) probe present in all models
    when TF protein translation model is used. Set to NULL (default)
    when translation model is not used.}
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
  \item{TFs}{The transcription factors that are tested by
    including them in the models one at a time.}
  \item{targets}{The target genes present in all models.}
}

\details{The models are formed by calling \code{\link{GPLearn}}.
  If there is no value given to the transcription
  factor, a model without protein translation is used.
  Without protein translation model, some
  known targets are needed. If known
  targets are given, a model is first created with only the transcription
  factor and the known targets. The parameters extracted from this model
  are used as initial parameters of the models with test targets.
  
  \code{GPRankTFs} is very similar to \code{GPRankTargets}, except
  it loops over candidate regulators, not candidate targets.
}

\value{
  The function returns a scoreList
  containing the genes, parameters and log-likelihoods of the models
  If returnModels is true,
  the function returns a list of the calculated models.
}
\author{Antti Honkela, Jonatan Ropponen, Magnus Rattray, Neil D. Lawrence}
\seealso{
  \code{\link{GPLearn}, \linkS4class{scoreList}, \link{generateModels}}.
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
}}
\keyword{model}
