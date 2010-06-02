gpOptions <- function (approx="ftc") {

  options = list()

  options$approx = approx
  
  ## Select type of optimiser.
  options$optimiser = "SCG"

  ## Set to true to learn output scales.
  options$learnScales = FALSE

  ## Set to true to scale outputs to variance 1.
  options$scale2var1 = FALSE

  ## Set to true to optimise beta.
  options$optimiseBeta = FALSE
  if (approx != "ftc")
    options$optimiseBeta = TRUE

  ## Set to a given mean function to have a mean function.
  options$meanFunction = list()
  ## Options structure for mean function options.
  options$meanFunctionOptions = list()

  ## Set to 1 if output processes have a shared variance.
  options$isSpherical = TRUE

  ## Set to 1 if there is data missing in the target matrix.
  options$isMissingData = FALSE

  if (options$approx == "ftc") {
    ## bog-standard kernel.
    ## The R version of the kern field is a bit more structured than the MATLAB version.
    options$kern = list(type="cmpnd",comp=list("rbf", "bias", "white"))
    options$numActive = 0
    options$beta = list()
  } else if (options$approx %in% c('fitc', 'pitc', 'dtc', 'dtcvar')) {
    options$kern = list(type="cmpnd",comp=list("rbf", "bias", "white"))
    options$numActive = 100
    options$beta = 1e+3
    ## Option to fix the inducing variables to other latent points.
    options$fixInducing = 0
    options$fixIndices = list()
  }

  options$computeS = FALSE
  return (options)
}