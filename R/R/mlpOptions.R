mlpOptions <- function(numHidden=20) {

  options$hiddenDim = numHidden
  options$activeFunc = 'linear'
  options$optimiser = 'SCG'

  return (options)
}
