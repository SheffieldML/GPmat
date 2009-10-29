\name{GPrankTargets}
\Rdversion{1.0}
\alias{GPrankTargets}
\title{Ranking possible target genes by forming optimized models}
\description{
  Ranks possible target genes by forming optimized models with a fixed transcription factor, a set of known target genes and targets to be tested. The transcription factor and the known targets are always included in the models while the tested targets are tested by including them in the models one at a time. The function determines itself whether to use GPsim or GPdisim based on the input arguments. Finally, the formed models are sorted based on the logarithmic value of their marginal likelihood.
}
\usage{
data <- GPrankTargets(preprocData, TF = NULL, knownTargets = NULL, testTargets = NULL, filterLimit = 1.5, useMedians = TRUE, returnScoreList = TRUE, returnModels = FALSE)
}
\arguments{
  \item{preprocData}{The preprocessed data to be used.}
  \item{TF}{The transcription factor present in all models.}
  \item{knownTargets}{The target genes present in all models.}
  \item{testTargets}{Target genes that are tested by including them in the models one at a time.}
  \item{filterLimit}{Genes with a filtering value above this figure are accepted after filtering. If this value is 0, all genes will be accepted.}
  \item{useMedians}{A logical value determining whether medians are used for filtering the genes; if false, means are used.}
  \item{returnScoreList}{A logical value determining whether the function returns a scoreList.}
  \item{returnModels}{A logical value determining whether the function returns the calculated models.}
}
\details{
  The models are formed by calling the function formModel, which in turn calls GPrank. If there is no value given to the transcription factor, GPsim is used. Otherwise, GPdisim is used. If GPsim is used, not giving known transcription factors results in an error. If known targets exist, a model is first created with only the transcription factor and the known targets. The parameters extracted from this model are used as initial parameters of the models with test targets.
}
\value{
  If returnScoreList is true, the function returns a scoreList containing the genes, parameters and log-likelihoods of the models as well as the value of useGPsim in each model. If returnModels is true, the function returns a list of the calculated models. If both are true, it returns both. The models and other values have been sorted in the order of the log-likelihoods, starting from the highest.
}
\author{
  Jonatan Ropponen
}
\seealso{
  \code{\link{GPrank}, \link{formModel}, \link{GPrankTFs}, \link{scoreList}, \link{generateModels}}.
}
\examples{
## missing
}
\keyword{model}
