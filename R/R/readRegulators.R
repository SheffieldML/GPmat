readRegulators <- function(fileName) {

  regulatorTable <- read.table(fileName)
  amountOfRegulators <- dim(regulatorTable)[1]
  regulators <- array(dim = amountOfRegulators)

  for (i in 1:amountOfRegulators) {
    currentRegulator <- as.character(regulatorTable[i, 1])
    tryCatch({
      regulators[i] <- as.character(mgu74av2ALIAS2PROBE[currentRegulator])
    }, error = function(err) {
      text = c("Value for gene", currentRegulator, "not found. \n")
      cat(text)
    })
  }

  return(regulators)
}


readPMIDRegulators <- function(fileName) {

  regulatorTable <- read.table(fileName)
  amountOfRegulators <- dim(regulatorTable)[1]
  regulators <- array()

  # counter for the amount of regulators found so far
  k <- 0

  for (i in 1:amountOfRegulators) {
    currentRegulator <- as.character(regulatorTable[i, 1])
    tryCatch({
      currentProbes <- as.character(mgu74av2PMID2PROBE[currentRegulator])
    }, error = function(err) {
      text = c("Value for gene", currentRegulator, "not found. \n")
      cat(text)
    })
    for (j in 1:length(currentProbes)) {
      regulators[k] <- currentProbes[j]
      k <- k + 1
    }
  }

  return(regulators)
}


readProbeRegulators <- function(fileName) {

  regulatorTable <- read.table(fileName)
  amountOfRegulators <- dim(regulatorTable)[1]
  regulators <- array(dim = amountOfRegulators)

  for (i in 1:amountOfRegulators) {
    regulators[i] <- as.character(regulatorTable[i, 1])
  }

  return(regulators)
}


removeDuplicates <- function(regulators) {

  for (i in 1:(length(regulators)-1)) {
    currentGene <- regulators[i]

    if (!is.na(currentGene)) {
      for (j in (i+1):length(regulators)) {
        if (!is.na(regulators[j]) && !is.na(charmatch(currentGene, regulators[j]))) {
	  regulators[j] <- NA
        }
      }
    }
  }

  # removing the NAs

  # new array of regulators
  regulators2 <- array()

  # counter of the genes in regulator2
  j <- 1

  for (i in 1:length(regulators)) {
    if (!is.na(regulators[i])) {
      regulators2[j] <- regulators[i]
      j <- j + 1
    }
  }

  return(regulators)
}
