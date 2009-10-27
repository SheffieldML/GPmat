generateModels <- function(preprocData = NULL, y = NULL, yvar = NULL, allGenes = NULL, times = NULL, scoreList = NULL, params = NULL, genes = NULL, useGPsim = NULL) {

  options(error = recover)
  options <- list(includeNoise=0, optimiser="CG")
  options$fix$index <- 4
  options$fix$value <- expTransform(c(1, 1), "xtoa")

  if (is.scoreList(scoreList)) {
    params <- scoreList@params
    genes <- scoreList@genes
    useGPsim <- scoreList@useGPsim
  }

  if (!is.null(preprocData)) {
    y <- preprocData$y
    yvar <- preprocData$yvar
    allGenes <- preprocData$genes
    times <- preprocData$times
  }

  listOfModels <- list()

  # for each model in scoreList
  for (i in 1:length(genes)) {

    currentGenes <- genes[[i]]
    currentParams <- params[[i]]
    useGPsimCurrently <- useGPsim[[i]]

    data <- searchExpressionData(y, yvar, allGenes, times, currentGenes)
    currentY <- data$y
    currentYvar <- data$yvar
    currentGenes <- data$genes

    #model <- gpsimCreate(Ngenes = length(currentGenes), Ntf = 1, times, currentY, currentYvar, options)
    #model <- gpsimExpandParam(model, currentParams)

    if (useGPsimCurrently) {
      model <- gpsimCreate(Ngenes = length(currentGenes), Ntf = 1, times, currentY, currentYvar, options)
      model <- gpsimExpandParam(model, currentParams)
    }
    else {
      model <- gpdisimCreate(Ngenes = length(currentGenes) - 1, Ntf = 1, times, currentY, currentYvar, options)
      model <- gpdisimExpandParam(model, currentParams)
    }

    listOfModels[[i]] <- model
  }

  return (listOfModels)
}


searchExpressionData <- function(y, yvar, allGenes, times, searchedGenes) {

  amountOfRows <- length(times)
  Nrep <- length(y)

  #counting the amount of found genes

  #counter for found genes
  k <- 0

  # for each specified gene
  for (i in 1:length(searchedGenes)) {
    searchedGene <- searchedGenes[i]

    # for each gene of the data
    for (j in 1:length(allGenes)) {
      if (!is.na(charmatch(searchedGene, allGenes[j]))) {
        k <- k + 1
      }
    }
  }

  foundGenes <- array(dim = k)
  foundY <- list()
  foundYvar <- list()

  foundY <- array(dim = c(amountOfRows, k))
  foundYvar <- array(dim = c(amountOfRows, k))

  #resetting the counter
  k <- 0

  #copying the data related to the found genes

  # for each specified gene
  for (i in 1:length(searchedGenes)) {
    searchedGene <- searchedGenes[i]

    # for each gene of the data
    for (j in 1:length(allGenes)) {
      if (!is.na(charmatch(searchedGene, allGenes[j]))) {
        k <- k + 1
        foundGenes[k] <- searchedGene
        for (l in 1: amountOfRows) {
          foundY[l, k] <- y[l, j]
          foundYvar[l, k] <- yvar[l, j]
	}
      }
    }
  }

  newData <- list(y = foundY, yvar = foundYvar, genes = foundGenes, times = times)
  return (newData)
}
