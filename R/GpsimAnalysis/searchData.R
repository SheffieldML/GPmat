searchData <- function(yFull, genes, pcts, searchedGenes, amountOfColumns) {

  #counting the amount of found genes

  #counter for found genes
  k <- 0

  # for each target gene
  for (i in 1:length(searchedGenes)) {
    searchedGene <- Genes[i]

    # for each gene of the data
    for (j in 1:length(genes)) {
      if (!is.na(charmatch(searchedGene, genes[j]))) {
        k <- k + 1
      }
    }
  }

  foundGenes <- array(dim = k)
  foundY <- array(dim = c(k, amountOfColumns))
  foundPcts <- array(dim = c(k, amountOfColumns, 5))

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
        for (l in 1: amountOfColumns) {
          foundY[k, l] <- yFull[j, l]
          # for each percentile
          for (m in 1:5) {
            foundPcts[k, l, m] <- pcts[j, l, m]
          }
        }
      }
    }
  }
  newData <- list(yFull = foundY, genes = foundGenes, pcts = foundPcts)
  return (newData)
}
