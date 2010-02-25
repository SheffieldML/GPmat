\name{formModel}
\Rdversion{1.0}
\alias{formModel}
\title{Forming optimized models with error handling.}
\description{
  Forms an optimized model by calling the function GPLearn with the given parameters. If GPLearn produces an error, it is called again several times with randomized parameters. The parameters and the log-likelihood of the model are also determined.
}
\usage{
returnData <- formModel(preprocData, TF = NULL, knownTargets = NULL, testTarget = NULL, useGPsim = TRUE, fixedParams = FALSE, initParams = NULL, fixComps = 1)
}
\arguments{
  \item{preprocData}{The preprocessed data to be used.}
  \item{TF}{The transcription factor of the model. In GPsim, there can be only one transcription factor for each model.}
  \item{knownTargets}{The known target genes of the model.}
  \item{testTarget}{The possible test target added to targets.}
  \item{useGPsim}{A logical value determining whether GPsim is used; if false, GPdisim is used.}
  \item{randomize}{A logical value determining whether the parameters of the model are randomized before optimization.}
  \item{addPriors}{A logical value determining whether priors are added to the model.}
  \item{search}{A logical value determining whether the preprocessed data is searched for the data of the transcription factor and target genes.}
  \item{fixedParams}{A logical value determining whether the initial parameters are fixed.}
  \item{initParams}{The initial parameters for the model. They are only used if fixedParams is TRUE.}
  \item{fixComps}{The blocks of the kern cache that are to be fixed.}
}
\value{
  Returns a list containing three items: the model returned by GPLearn, the parameters of the model and the log-likelihood of the model.
}
\author{
  Jonatan Ropponen
}
\seealso{
\code{\link{GPLearn}, \link{GPrankTargets}, \link{GPrankTFs}}.
}
\examples{
## missing
}
\keyword{model}
