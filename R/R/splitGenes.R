splitGenesFixedTF <- function(genes, TF = "100001_at", splits = 1) {

  amountOfGenes <- length(genes)

  geneList <- list()

  for (i in 1:splits) {

    searchedGenes <- array()
    searchedGenes[1] <- TF

    # counter of genes in searchedGenes
    k <- 2

    for (j in (floor(amountOfGenes/splits*(i-1))+1):floor(amountOfGenes/splits*i)) {
      if (!is.na(genes[j]) && is.na(charmatch(TF, genes[j]))) {
        searchedGenes[k] <- genes[j]
        k <- k + 1
      }
    }
    geneList[[i]] <- searchedGenes
  }

  return(geneList)
}


splitGenesFixedTargets <- function(genes, targets = "100001_at", splits = 1) {

  amountOfGenes <- length(genes)

  geneList <- list()

  for (i in 1:splits) {

    searchedGenes <- array()
    for (i in 1:length(targets)) {
      searchedGenes[i] <- targets[i]
    }

    # counter of genes in searchedGenes
    k <- length(targets) + 1

    for (j in (floor(amountOfGenes/splits*(i-1))+1):floor(amountOfGenes/splits*i)) {
      if (!is.na(genes[j]) && is.na(charmatch(TF, genes[j]))) {
        searchedGenes[k] <- genes[j]
        k <- k + 1
      }
    }
    geneList[[i]] <- searchedGenes
  }

  return(geneList)
}

writeGeneList <- function(geneList, fileName = "genes") {

  for (i in 1:length(geneList)) {
    genes <- geneList[[i]]
    fileName2 <- paste(fileName, i, sep="")
      text <- ""
    for (i in 1:length(genes)) {
      text <- c(text, genes[i])
    }
    write(text, file = fileName2)
  }
}
