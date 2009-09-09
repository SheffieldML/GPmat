\name{GPscoreListFixedTF}
\Rdversion{1.0}
\alias{GPscoreListFixedTF}
\title{Forming optimized models with a fixed transcription factor}
\description{
  Forms optimized models with a fixed transcription factor, a set of known target genes and targets to be tested. The transcription factor and the known targets are always included in the models while the tested targets are tested by including them in the models one at a time. The function determines itself whether to use GPsim or GPdisim based on the input arguments. Finally, the formed models are sorted based on the logarithmic value of their marginal likelihood.
}
\usage{
scoreList <- GPscoreListFixedTF(preprocData, TF = NULL, knownTargets = NULL, testTargets = NULL)
}
\arguments{
  \item{preprocData}{The preprocessed data to be used.}
  \item{TF}{The transcription factor present in all models.}
  \item{knownTargets}{The target genes present in all models.}
  \item{testTargets}{Target genes that are tested by including them in the models one at a time.}
}
\details{
  The models are formed by calling the function formModel, which in turn calls GPrank. If there is no value given to the transcription factor, GPsim is used. Otherwise, GPdisim is used. If GPsim is used, not giving known transcription factors results in an error. If known targets exist, a model is first created with only the transcription factor and the known targets. The parameters extracted from this model are used as initial parameters of the models with test targets.
}
\value{
  Returns a list containing the data, parameters and log-likelihoods of the models. The values have been sorted in the order of the log-likelihoods, starting from the highest.
}
\seealso{
  \code{\link{GPrank}, \link{formModel}}.
}
\examples{
## missing
}
\keyword{model}
