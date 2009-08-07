readRegulators <- function(fileName) {

  regulatorTable <- read.table(fileName)
  amountOfRegulators <- dim(regulatorTable)[1]
  regulators <- array(dim = amountOfRegulators)

  for (i in 1:amountOfRegulators) {
    currentRegulator <- as.character(regulatorTable[i, 1])
    tryCatch({
      regulators[i] <- as.character(mgu74av2ALIAS2PROBE[currentRegulator])
    }, error = function(err) {
      text = c("Value for gene ", currentRegulator, " not found. \n")
      cat(text)
    })
  }

  return(regulators)
}
