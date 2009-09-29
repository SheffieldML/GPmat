\name{GPrankTFs}
\Rdversion{1.0}
\alias{GPrankTFs}
\title{Ranking possible transcription factors by forming optimized models}
\description{
  Ranks possible trancription factors by forming optimized models with fixed target genes and a set of transcription factors to be tested. The targets are always included in the models while the transcription factors are tested by including them in the models one at a time. The function uses GPsim when a transcription factor is not used. Finally, the formed models are sorted based on the logarithmic value of their marginal likelihood.
}
\usage{
scoreList <- GPrankTFs(preprocData, TFs = NULL, targets = NULL, filterLimit = 1.5, useMedians = TRUE)
}
\arguments{
  \item{preprocData}{The preprocessed data to be used.}
  \item{TFs}{The transcription factors that are tested by including them in the models one at a time.}
  \item{targets}{The target genes present in all models.}
  \item{filterLimit}{Genes with a filtering value above this figure are accepted after filtering. If this value is 0, all genes will be accepted.}
  \item{useMedians}{A logical value determining whether medians are used for filtering the genes; if false, means are used.}
}
\details{
  The models are formed by calling the function formModel, which in turn calls GPrank. If there is no value given to the transcription factor, GPsim is used. Otherwise, GPdisim is used. If GPsim is used, not giving known transcription factors results in an error. The initial model is created with only the target genes. The parameters extracted from this model are used as initial parameters of the models with test targets.
}
\value{
  Returns a list containing the data, parameters and log-likelihoods of the models. The values have been sorted in the order of the log-likelihoods, starting from the highest.
}
\seealso{
  \code{\link{GPrank}, \link{formModel}, \link{GPrankTargets}}.
}
\examples{
## missing
}
\keyword{model}
