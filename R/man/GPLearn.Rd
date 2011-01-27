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
  model <- GPLearn(preprocData, TF = NULL, targets = NULL,
    useGpdisim = !is.null(TF), randomize = FALSE, addPriors = FALSE,
    fixedParams = FALSE, initParams = NULL, initialZero = TRUE,
    fixComps = NULL, dontOptimise = FALSE,
    allowNegativeSensitivities = FALSE, quiet = FALSE,
    gpsimOptions = NULL, allArgs = NULL)
}
\arguments{
  \item{preprocData}{The preprocessed data to be used.}
  \item{TF}{The probe corresponding to the transcription factor (TF) mRNA
    if TF protein translation model is used, or NULL (default) if the
    translation model is not used.}
  \item{targets}{The target genes of the model.}
  \item{useGpdisim}{A logical value determining whether a model of
    translation is included.  By default TRUE if TF is set, FALSE
    if TF is unset.}
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
  \item{allowNegativeSensitivities}{Allow sensitivities to go negative.
    This is an experimental feature, and the negative values have no
    physical interpretation.}
  \item{quiet}{Suppress optimiser output.}
  \item{gpsimOptions}{Internal: additional options to pass to
    gp[di]simCreate.}
  \item{allArgs}{A list of arguments that can be used to override
    ones with the same name.}
}
\value{
  Returns the optimized model.
}
\author{Antti Honkela, Pei Gao, Jonatan Ropponen, Magnus Rattray, Neil D. Lawrence}
\seealso{
  \code{\link{GPRankTargets}, \link{GPRankTFs}}.
}
\examples{
  # Load a mmgmos preprocessed fragment of the Drosophila developmental
  # time series
  data(drosophila_gpsim_fragment)

  # Get the target probe names
  library(annotate)
  aliasMapping <- getAnnMap("ALIAS2PROBE",
                    annotation(drosophila_gpsim_fragment))
  twi <- get('twi', env=aliasMapping)
  fbgnMapping <- getAnnMap("FLYBASE2PROBE",
                   annotation(drosophila_gpsim_fragment))
  targetProbe <- get('FBgn0035257', env=fbgnMapping)

  # Create the model but do not optimise (rarely needed...)
  model <- GPLearn(drosophila_gpsim_fragment,
                   TF=twi, targets=targetProbe,
                   useGpdisim=TRUE, quiet=TRUE,
                   dontOptimise=TRUE)
\dontrun{
  # Create and learn the model
  model <- GPLearn(drosophila_gpsim_fragment,
                   TF=twi, targets=targetProbe,
                   useGpdisim=TRUE, quiet=TRUE)
}}
\keyword{model}
