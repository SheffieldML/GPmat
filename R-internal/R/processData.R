processData <- function(data, replicates = NULL, searchGenes = NULL) {

  require(puma)

  yFull <- exprs(data)
  times <- pData(data)
  genes <- rownames(exprs(data))

  times <- times[,grep('time', colnames(times))]

  numberOfRows <- length(genes)
  numberOfColumns <- length(colnames(exprs(data)))

  if (is.null(replicates))
    replicates <- rep(1, numberOfColumns)
  
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

  # Search specified genes if requested
  if (! is.null(searchGenes) ) {
    newData <- searchData(yFull, genes, pcts, searchGenes, numberOfColumns)
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
  rownames(yFullVar) <- rownames(yFull)
  colnames(yFullVar) <- colnames(yFull)
  
  y <- list()
  yvar <- list()

  repids <- unique(replicates)
  for (i in seq(along=repids)) {
    y[[i]] <- yFull[replicates==repids[i],]
    yvar[[i]] <- yFullVar[replicates==repids[i],]
  }
  times <- as.array(times[replicates==repids[1]])
  ##times <- array(times[, 1], dim = c(numberOfColumns)) 

  zScores <- as.array(colMeans(yFull / sqrt(yFullVar)))

  return (GPdata(y = y, yvar = yvar, times = times, genes = genes,
                 scale = as.array(scale), zScores = zScores,
                 experiments=array(), annotation=annotation(data),
                 phenoData=phenoData(data), featureData=featureData(data)))
}


processStructuredData <- function(data, experiments, searchGenes = NULL) {

  require(puma)

  yFull <- exprs(data)
  times <- pData(data)
  genes <- rownames(exprs(data))

  times <- times[,grep('time', colnames(times))]

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

  # Search specified genes if requested
  if (! is.null(searchGenes) ) {
    newData <- searchData(yFull, genes, pcts, searchGenes, numberOfColumns)
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
  rownames(yFullVar) <- rownames(yFull)
  colnames(yFullVar) <- colnames(yFull)
  
  y <- list()
  yvar <- list()

  y[[1]] <- yFull
  yvar[[1]] <- yFullVar
  times <- as.array(times)

  zScores <- as.array(colMeans(yFull / sqrt(yFullVar)))

  return (GPdata(y = y, yvar = yvar, times = times, genes = genes,
                 scale = as.array(scale), zScores = zScores,
                 experiments=as.array(experiments),
                 annotation=annotation(data), phenoData=phenoData(data),
                 featureData=featureData(data)))
}


processRawData <- function(data, times, replicates=NULL) {
  yFull <- as.matrix(data)
  genes <- rownames(data)

  numberOfRows <- length(genes)
  numberOfColumns <- length(colnames(data))

  if (is.null(replicates))
    replicates <- rep(1, numberOfColumns)
  
  normalisation <- colMeans(yFull) - mean(yFull)
  zeroArray <- array(0, dim=c(numberOfRows, numberOfColumns))

  # The default operation of sweep is "-".
  yFull <- yFull + sweep(zeroArray, 2, normalisation)

  yFull <- t(exp(yFull))

  scale <- sd(yFull)
  scaleMat <- array(1, dim = c(numberOfColumns, 1)) %*% scale
  yFull <- yFull / scaleMat

  y <- list()
  yvar <- list()

  repids <- unique(replicates)
  for (i in seq(along=repids)) {
    y[[i]] <- yFull[replicates==repids[i],]
    yvar[[i]] <- 0 * y[[i]]
  }
  times <- as.array(times[replicates==repids[1]])

  zScores <- array(NA, dim=length(genes))
  names(zScores) <- colnames(y[[1]])
  
  return (GPdata(y = y, yvar = yvar, times = times, genes = genes,
                 scale = as.array(scale), zScores = zScores,
                 experiments=array(), annotation="NA",
                 phenoData=new("AnnotatedDataFrame"),
                 featureData=new("AnnotatedDataFrame")))
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
