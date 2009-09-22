processData <- function(data, searchedGenes = "100001_at", search = FALSE) {

  require(puma)

  yFull <- exprs(data)
  times <- pData(data)
  genes <- rownames(exprs(data))

  amountOfRows <- length(genes)
  amountOfColumns <- length(colnames(exprs(data)))

  # Normalisation is done for all genes, before specific genes are taken from the array.
  normalisation <- array(0, dim = c(1, amountOfColumns))

  for (i in 1:amountOfColumns) {
    normalisation[i] <- mean(yFull[, i]) - mean(yFull)
  }
  zeroArray <- array(0, dim=c(amountOfRows, amountOfColumns))

  # The default operation of sweep is "-".
  yFull <- yFull + sweep(zeroArray, 2, normalisation)

  pcts <- array(dim = c(amountOfRows, amountOfColumns, 5))

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
    newData <- searchData(yFull, genes, pcts, searchedGenes, amountOfColumns)
    yFull <- newData$yFull
    genes <- newData$genes
    pcts <- newData$pcts
  }

  amountOfRows <- length(genes)
  amountOfColumns <- length(colnames(exprs(data)))

  yFullVar <- array(dim = c(amountOfRows, amountOfColumns))

  for (i in 1:amountOfRows) {
    prof <- pcts[i,,]
    for (j in 1:amountOfColumns) {
      t <- distfit(exp(prof[j, ]), 'normal')
      yFull[i, j] <- t$par[1]
      yFullVar[i, j] <- t$par[2] ^ 2
    }
  }

  yFull <- t(yFull)
  yFullVar <- t(yFullVar)

  # Rescale so that average standard deviation of curves is 1.
  scale <- sd(yFull)
  scaleMat <- array(1, dim = c(amountOfColumns, 1)) %*% scale
  yFull <- yFull / scaleMat
  yFullVar <- yFullVar / scaleMat^2

  y <- list()
  y[[1]] <- yFull
  yvar <- list()
  yvar[[1]] <- yFullVar
  times <- array(times[, 1], dim = c(amountOfColumns)) 

  ratioData <- averageToSDRatio(usePreprocData = FALSE, y = y, yvar = yvar, times = times, genes = genes)

  preprocData <- list(y = y, yvar = yvar, times = times, genes = genes, scale = scale, ratioData = ratioData)

  return(preprocData)
}


searchData <- function(yFull, genes, pcts, searchedGenes, amountOfColumns) {

  #counting the amount of found genes

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
