processData <- function(data, searchedGenes = "100001_at", search = FALSE) {

  require(puma)

  yFull <- exprs(data)
  times <- pData(data)
  genes <- rownames(exprs(data))

  numberOfRows <- length(genes)
  numberOfColumns <- length(colnames(exprs(data)))

  # Normalisation is done for all genes, before specific genes are taken from the array.
  normalisation <- array(0, dim = c(1, numberOfColumns))

  for (i in 1:numberOfColumns) {
    normalisation[i] <- mean(yFull[, i]) - mean(yFull)
  }
  zeroArray <- array(0, dim=c(numberOfRows, numberOfColumns))

  # The default operation of sweep is "-".
  yFull <- yFull + sweep(zeroArray, 2, normalisation)

  pcts <- array(dim = c(numberOfRows, numberOfColumns, 5))

  # Loading the percentiles 5, 25, 50, 75 and 95.
  pcts[,,1] <- prcfive(data)
  pcts[,,2] <- prctwfive(data)
  pcts[,,3] <- prcfifty(data)
  pcts[,,4] <- prcsevfive(data)
  pcts[,,5] <- prcninfive(data)

  # normalising the percentiles
  for (i in 1:5) {
    pcts[,,i] <- pcts[,,i] + sweep(zeroArray, 2, normalisation)
  }

  # The variable search determines whether the specified genes are searched.
  if (search) {
    newData <- searchData(yFull, genes, pcts, searchedGenes, numberOfColumns)
    yFull <- newData$yFull
    genes <- newData$genes
    pcts <- newData$pcts
  }

  numberOfRows <- length(genes)
  numberOfColumns <- length(colnames(exprs(data)))

  yFullVar <- array(dim = c(numberOfRows, numberOfColumns))

  for (i in 1:numberOfRows) {
    prof <- pcts[i,,]
    for (j in 1:numberOfColumns) {
      t <- distfit(exp(prof[j, ]), 'normal')
      yFull[i, j] <- t$par[1]
      yFullVar[i, j] <- t$par[2] ^ 2
    }
  }

  yFull <- t(yFull)
  yFullVar <- t(yFullVar)

  # Rescale so that average standard deviation of curves is 1.
  scale <- sd(yFull)
  scaleMat <- array(1, dim = c(numberOfColumns, 1)) %*% scale
  yFull <- yFull / scaleMat
  yFullVar <- yFullVar / scaleMat^2

  y <- list()
  y[[1]] <- yFull
  yvar <- list()
  yvar[[1]] <- yFullVar
  times <- array(times[, 1], dim = c(numberOfColumns)) 

  ratioData <- averageToSDRatio(usePreprocData = FALSE, y = y, yvar = yvar, times = times, genes = genes)

  preprocData <- list(y = y, yvar = yvar, times = times, genes = genes, scale = scale, ratioData = ratioData)

  return(preprocData)
}


searchData <- function(yFull, genes, pcts, searchedGenes, numberOfColumns) {

  #counting the number of found genes

  #counter for found genes
  k <- 0

  # for each target gene
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
  foundY <- array(dim = c(k, numberOfColumns))
  foundPcts <- array(dim = c(k, numberOfColumns, 5))

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
        for (l in 1: numberOfColumns) {
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


filterGenes <- function(ratioData, filterLimit = 1.5, useMedians = TRUE) {

  genes <- ratioData$genes

  approvedGenes <- array()

  # counter for the next approved gene
  j <- 1

  for (i in 1:length(genes)) {
    if (useMedians) {
      if (ratioData$medians[[1]][i] >= filterLimit) {
	approvedGenes[j] <- genes[i]
	j <- j + 1
      }
    }
    else {
      if (ratioData$means[[1]][i] >= filterLimit) {
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

  numberOfRows <- length(times)
  Nrep <- length(y)

  ratios <- list()

  # for each gene of the data
  for (m in 1:Nrep) {
    ratios[[m]] <- array(dim = dim(y[[m]]))
    for (j in 1:length(genes)) {
      for (i in 1: numberOfRows) {
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
