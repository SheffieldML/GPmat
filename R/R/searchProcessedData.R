searchProcessedData <- function(preprocData, searchedGenes) {

  y <- preprocData$y
  yvar <- preprocData$yvar
  times <- preprocData$times
  genes <- preprocData$genes
  scale <- preprocData$scale

  amountOfRows <- length(times)
  Nrep <- length(y)

  #counting the amount of found genes

  #counter for found genes
  k <- 0

  # for each specified gene
  for (i in 1:length(searchedGenes)) {
    searchedGene <- searchedGenes[i]

    # for each gene of the data
    for (j in 1:length(genes)) {
      if (!is.na(charmatch(searchedGene, genes[j]))) {
        k <- k + 1
      }
    }
  }

  foundGenes <- array(dim = k)
  foundY <- list()
  foundYvar <- list()
  for (m in 1:Nrep) {
    foundY[[m]] <- array(dim = c(amountOfRows, k))
    foundYvar[[m]] <- array(dim = c(amountOfRows, k))
  }

  #resetting the counter
  k <- 0

  #copying the data related to the found genes

  # for each specified gene
  for (i in 1:length(searchedGenes)) {
    searchedGene <- searchedGenes[i]

    # for each gene of the data
    for (j in 1:length(genes)) {
      if (!is.na(charmatch(searchedGene, genes[j]))) {
        k <- k + 1
        foundGenes[k] <- searchedGene
	for (m in 1:Nrep) {
          for (l in 1: amountOfRows) {
            foundY[[m]][l, k] <- y[[m]][l, j]
            foundYvar[[m]][l, k] <- yvar[[m]][l, j]
          }
	}
      }
    }
  }
  newData <- list(y = foundY, yvar = foundYvar, genes = foundGenes, times = times, scale = scale)
  return (newData)
}
