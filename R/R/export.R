export.scores <- function(scores, datasetName='', experimentSet='',
    databaseFile='database.sqlite', preprocData=NULL, models=NULL,
    figpath=NULL, aliasTypes=c("SYMBOL", "GENENAME", "ENTREZID"),
    datasetSource='', datasetDescription='',
    datasetSaveLocation='', datasetFigureFilename='',
    experimentTimestamp=as.character(Sys.Date()),
    figureDesc='', figurePrio=0, regulator=NULL) {

  if (datasetName == '') {
    tryCatch(datasetName <- datasetName(scores), error=function(e) message('old version of scoreList object detected'))
  }

  if (experimentSet == '') {
    tryCatch(experimentSet <- experimentSet(scores), error=function(e) message('old version of scoreList object detected'))
  }

  if (datasetName == '') {
    stop('datasetName not defined')
  }

  if (experimentSet == '') {
    stop('experimentSet not defined')
  }

  experimentProducer <- paste('tigre', sessionInfo('tigre')[["otherPkgs"]][['tigre']]$Version)

  if (TF(scores) == "(see genes)") { # TODO
    loopTarget <- modelArgs(scores)[[1]]$targets
    loopVariable <- 1 # loop TFs, target set constant
    TF <- NULL
    targetSet <- modelArgs(scores)[[1]]$targets
  } else {
    loopTarget <- TF(scores)
    loopVariable <- 2 # loop targets, TF constant
    if (! is.null(regulator)) {
      TF <- regulator
    } else {
      TF <- modelArgs(scores)[[1]]$TF
      if (is.null(TF))
        TF <- "N/A"
    }
    targetSet <- NULL
  }

  db <- .openConnection(databaseFile)
  .createTables(db)

  loopVariables <- .getLoopVariables(db)
  if (length(loopVariables) > 0 && loopVariables[[1]] != loopVariable) {
    stop('Error: Same database cannot contain results from looping over TFs and looping over targets.')
  }

  # Generate models if necessary and possible
  #if (is.null(models) && !is.null(preprocData)) {
  #  message("Generating models...")
  #  models <- generateModels(preprocData, scores)
  #}

  figureData <- list()
  if (!is.null(figpath) && (!is.null(models) || !is.null(preprocData))) {
    message("Generating figures...")
    figureFilename <- file.path(figpath, paste('model_', loopTarget, '_${probe_name}.png', sep=""))
    len <- length(scores)
    for (k in seq(along=scores)) {
      message(paste("figure", k, "/", len))
      if (is.null(models))
        m <- generateModels(preprocData, scores[k])[[1]]
      else
        m <- models[[k]]
      png(file.path(figpath, paste('model_', loopTarget, '_', genes(scores)[k], '.png', sep="")))
      GPPlot(m)
      dev.off()
    }
  } else if (!is.null(models) || !is.null(preprocData)) {
    message("Generating figures...")
    figureFilename <- ''
    len <- length(scores)
    for (k in seq(along=scores)) {
      message(paste("figure", k, "/", len))
      if (is.null(models))
        m <- generateModels(preprocData, scores[k])[[1]]
      else
        m <- models[[k]]
      filename <- tempfile()
      png(filename)
      GPPlot(m)
      dev.off()
      figureData[[genes(scores)[[k]]]] <- filename
    }
  }

  if (modelArgs(scores)[[1]]$useGpdisim) {
    modelTranslation <- TRUE
    experimentName <- "GPDISIM"
  } else {
    modelTranslation <- FALSE
    experimentName <- "GPSIM"
  }

  experimentDesc <- experimentName
  if (length(knownTargets(scores)) > 1 || nchar(knownTargets(scores)) > 0) {
    experimentName <- paste(paste(experimentName, ',', sep=''), length(knownTargets(scores)), "known target(s)")
  }

  has.annotation <- FALSE
  if (!is.null(preprocData)) {
    has.annotation <- length(annotation(preprocData)) > 0
  }
  if (!has.annotation) {
    message('Warning: annotation(preprocData) is not defined, cannot add aliases')
  }

  datasetSpecies <- ''
  datasetPlatform <- ''
  if (!is.null(preprocData) && has.annotation) {
    datasetSpecies <- getAnnMap('ORGANISM', annotation(preprocData))
    datasetPlatform <- annotation(preprocData)
  }

  datasetId <- .addAndGetDatasetId(db, datasetName, datasetSpecies, datasetSource, datasetPlatform, datasetDescription, datasetSaveLocation, datasetFigureFilename)

  # Add aliases
  if (!is.null(preprocData) && has.annotation) {
    message("Adding aliases...")
    if (TF %in% featureNames(preprocData))
      probes <- c(names(loglikelihoods(scores)), TF, targetSet)
    else
      probes <- c(names(loglikelihoods(scores)), targetSet)
    len <- length(aliasTypes)
    for (i in seq_along(aliasTypes)) {
      message(paste("alias type", i, "/", len, paste("(", aliasTypes[[i]],")", sep="")))
      .addAliases(db, probes, preprocData, datasetId, aliasTypes[[i]])
    }
  }

  # Add z-scores
  if (!is.null(preprocData) && 'var.exprs' %in% assayDataElementNames(preprocData)) {
    message("Adding z-scores...")
    if (TF %in% featureNames(preprocData))
      probes <- c(names(loglikelihoods(scores)), TF, targetSet)
    else
      probes <- c(names(loglikelihoods(scores)), targetSet)
    
    data <- preprocData[probes,]
    y <- assayDataElement(data, 'exprs')
    yvar <- assayDataElement(data, 'var.exprs')

    zScores <- rowMeans(y / sqrt(yvar))
    .addZscores(db, datasetId, zScores)
  }

  if (!is.null(TF)) {
    regId <- .addAndGetRegulatorId(db, TF, datasetId)
  } else {
    regId <- NA
  }

  parameters <- params(scores)
  numberOfParameters <- length(parameters[[1]])
  parameterNames <- NA

  if (!is.null(models)) {
    parameterNames <- paste(names(modelExtractParam(models[[1]], only.values=FALSE)), collapse=',')
  } else if (!is.null(preprocData)) {
    m <- generateModels(preprocData, scores[1])
    parameterNames <- paste(names(modelExtractParam(m[[1]], only.values=FALSE)), collapse=',')
  }

  message("Adding results...")
  experimentId <- .addAndGetExperimentId(db, experimentName, experimentDesc, datasetId, regulatorId=regId, loopVariable=loopVariable, modelTranslation=modelTranslation, numberOfParameters=numberOfParameters, parameterNames=parameterNames, producer=experimentProducer, timestamp=experimentTimestamp)
  .addResults(db, experimentId, loglikelihoods(scores), baseloglikelihoods(scores), parameters)

  if (!is.null(figpath) || length(figureData) > 0) {
    message("Adding figures...")
    .addFigures(db, experimentId, figureFilename, name=paste(experimentName, loopTarget), description=figureDesc, priority=figurePrio, figureData=figureData)
  }

  if (!is.null(targetSet)) {
    .addTargetSet(db, experimentId, targetSet)
  }

  rootId <- .addAndGetExperimentSetId(db, 'All experiments', NA)
  setId <- .addAndGetExperimentSetId(db, paste(experimentSet, ' (', datasetName, ')', sep=''), rootId)
  .addExperimentSetExperiments(db, setId, experimentId)

  dbDisconnect(db)
  message("All done.")
}


