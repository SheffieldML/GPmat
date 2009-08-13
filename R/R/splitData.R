splitDataFixedTF <- function(preprocData, TF = "100001_at", splits = 1) {

  genes <- preprocData$genes
  amountOfGenes <- length(genes)

  newData <- list()
  searchedGenes <- array()
  searchedGenes[1] <- TF

  # counter of genes in searchedGenes
  k <- 2

  for (i in 1:splits) {
    
    for (j in (floor(amountOfGenes/splits*(i-1))+1):floor(amountOfGenes/splits*i)) {
      if (!is.na(genes[j]) && is.na(charmatch(TF, genes[j]))) {
        searchedGenes[k] <- genes[j]
        k <- k + 1
      }
    }

    k <- 2
    newData[[i]] <- searchProcessedData(preprocData, searchedGenes)
  }

  return(newData)
}


splitDataFixedTargets <- function(preprocData, targets = "100001_at", splits = 1) {

  genes <- preprocData$genes
  amountOfGenes <- length(genes)

  newData <- list()
  searchedGenes <- array()

  for (i in 1:length(targets)) {
    searchedGenes[i] <- targets[i]
  }

  # counter of genes in searchedGenes
  k <- length(targets) + 1

  for (i in 1:splits) {
    
    for (j in (floor(amountOfGenes/splits*(i-1))+1):floor(amountOfGenes/splits*i)) {
      if (!is.na(genes[j]) && is.na(charmatch(TF, genes[j]))) {
        searchedGenes[k] <- genes[j]
        k <- k + 1
      }
    }

    k <- 2
    newData[[i]] <- searchProcessedData(preprocData, searchedGenes)
  }

  return(newData)
}
