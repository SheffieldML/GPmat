\name{GPrank}
\Rdversion{1.0}
\alias{GPrank}
\title{Forming optimized models}
\description{
  Forms an optimized model of the desired genes. The function can form a model with GPsim or GPdisim and it's also possible to use initial parameters or fix parameters for future use. The genes can also be filtered based on ratios calculated from the expression values. The given data can also be searched for the data of specific genes.
}
\usage{
data <- GPrank(preprocData, TF = NULL, targets = NULL, useGPsim = FALSE, filterLimit = 0, useMedians = TRUE, randomize = FALSE, addPriors = FALSE, fixedParams = FALSE, initParams = NULL, fixComps = 1)
}
\arguments{
  \item{preprocData}{The preprocessed data to be used.}
  \item{TF}{The transcription factor of the model. In GPsim, there can be only one transcription factor for each model.}
  \item{targets}{The target genes of the model.}
  \item{useGPsim}{A logical value determining whether GPsim is used; if false, GPdisim is used.}
  \item{filterLimit}{Genes with a filtering value above this figure are accepted after filtering. If this value is 0, all genes will be accepted.}
  \item{useMedians}{A logical value determining whether medians are used for filtering the genes; if false, means are used.}
  \item{randomize}{A logical value determining whether the parameters of the model are randomized before optimization.}
  \item{addPriors}{A logical value determining whether priors are added to the model.}
  \item{search}{A logical value determining whether the preprocessed data is searched for the data of the transcription factor and target genes.}
  \item{fixedParams}{A logical value determining whether the initial parameters are fixed.}
  \item{initParams}{The initial parameters for the model. They are only used if fixedParams is TRUE.}
  \item{fixComps}{The blocks of the kern cache that are to be fixed.}
}
\value{
  Returns the data containing the model itself and the genes used.
}
\seealso{
  \code{\link{gpsimCreate}, \link{GPrankTargets}, \link{GPrankTFs}, \link{filterGenes}}.
}
\examples{
## missing
}
\keyword{model}
