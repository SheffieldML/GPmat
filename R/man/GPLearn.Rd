\name{GPLearn}
\Rdversion{1.0}
\alias{GPLearn}
\title{Fit a GP model}
\description{
  Forms an optimized model of the desired genes. The function can form a
  model with GPsim or GPdisim and it's also possible to use initial
  parameters or fix parameters for future use. The genes can also be
  filtered based on ratios calculated from the expression values. The
  given data can also be searched for the data of specific genes.
}
\usage{
GPLearn <- function(preprocData, TF = NULL, targets = NULL, useGpdisim = FALSE, randomize = FALSE, addPriors = FALSE, fixedParams = FALSE, initParams = NULL, initialZero = TRUE, fixComps = NULL, dontOptimise = FALSE, gpsimOptions = NULL, allArgs = NULL)
}
\arguments{
  \item{preprocData}{The preprocessed data to be used.}
  \item{TF}{The transcription factor of the model. In GPsim, there can
    be only one transcription factor for each model.}
  \item{targets}{The target genes of the model.}
  \item{useGpdisim}{A logical value determining whether GPDISIM is used,
    FALSE (default) implies GPSIM.}
  \item{randomize}{A logical value determining whether the parameters
    of the model are randomized before optimization.}
  \item{addPriors}{A logical value determining whether priors are
    added to the model.}
  \item{fixedParams}{A logical value determining whether the initial
    parameters are fixed.}
  \item{initParams}{The initial parameters for the model. In
    combination with fixedParams a value NA denotes parameters to
    learn.}
  \item{initialZero}{Assume a zero initial TF protein concentration,
    default = TRUE.}
  \item{fixComps}{The blocks of the kernel the parameters of which are
    to be fixed. To be used together with fixedParams and initParams.}
  \item{dontOptimise}{Just create the model, do not run optimisation.}
  \item{gpsimOptions}{Internal: additional options to pass to
    gp[di]simCreate.}
  \item{allArgs}{A list of arguments that can be used to override
    ones with the same name.}
}
\value{
  Returns the optimized model.
}
\author{
  Jonatan Ropponen, Antti Honkela
}
\seealso{
  \code{\link{gpsimCreate}, \link{GPRankTargets}, \link{GPRankTFs}}.
}
\examples{
## missing
}
\keyword{model}
