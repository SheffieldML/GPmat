filterGenes <- function(ratioData, filterLimit = 1.5, useMedians = TRUE) {

  genes <- ratioData$genes

  approvedGenes <- array()

  # counter for the next approved gene
  j <- 1

  for (i in 1:length(genes)) {
    if (useMedians) {
      if (ratioData$medians[[1]][i] > filterLimit) {
	approvedGenes[j] <- genes[i]
	j <- j + 1
      }
    }
    else {
      if (ratioData$means[[1]][i] > filterLimit) {
	approvedGenes[j] <- genes[i]
	j <- j + 1
      }
    }
  }

  return (approvedGenes)
}



averageToSDRatio <- function(y = NULL, yvar = NULL, times = NULL, genes = NULL, usePreprocData = FALSE, preprocData = NULL, searchedGenes = "100001_at", search = FALSE) {

  if (usePreprocData) {
    if (search) {
      newData <- searchProcessedData(preprocData, searchedGenes)
    }
    else {
      newData <- preprocData
    }

    y <- newData$y
    yvar <- newData$yvar
    times <- newData$times
    genes <- newData$genes
  }

  amountOfRows <- length(times)
  Nrep <- length(y)

  ratios <- list()

  # for each gene of the data
  for (m in 1:Nrep) {
    ratios[[m]] <- array(dim = dim(y[[m]]))
    for (j in 1:length(genes)) {
      for (i in 1: amountOfRows) {
        average <- y[[m]][i, j]
        variance <- yvar[[m]][i, j]
	ratios[[m]][i, j] <- average / sqrt(variance)
      }
    }
  }

  means <- list()
  medians <- list()

  # counting the mean and median of ratios
  for (m in 1:Nrep) {
    means[[m]] <- array(dim = length(genes))
    medians[[m]] <- array(dim = length(genes))
    for (j in 1:length(genes)) {
      means[[m]][j] <- mean(ratios[[m]][, j])
      medians[[m]][j] <- median(ratios[[m]][, j])
    }
  }

  ratioData <- list()
  ratioData$ratios <- ratios
  ratioData$means <- means
  ratioData$medians <- medians
  ratioData$genes <- genes

  return (ratioData)
}



writeRatioData <- function(ratioData, genes, fileName) {

  text <- ""
  for (i in 1:length(genes)) {
    gene <- genes[i]
    #ratios <- ratioData$ratios[[1]][, i]
    mean <- ratioData$means[[1]][i]
    median <- ratioData$medians[[1]][i]
    #newText <- c(gene, "\n", "ratios: ", ratios, "\n", 
#"mean: ", mean, "\n", "median: ", median, "\n")
    newText <- c(gene, "mean: ", mean, "median: ", median, "\n")
    text <- c(text, newText)
  }
  write(text, file = fileName)

}
