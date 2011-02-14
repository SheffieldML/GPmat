require(gplots)

GPPlot <- function(data, savepath = '', nameMapping = NULL,
                   predt = NULL, fileOutput=FALSE, plotTime=NULL) {
  FONTSIZE <- 10;
  LINEWIDTH <- 1;
  MARKERSIZE <- 10;

  if (is.list(data))
    data <- data[[1]]
  
  if (is.GPModel(data))
    model <- modelStruct(data)
  else {
    if ("model" %in% names(data))
      model <- data$model
    else
      model <- data
  }

  if ("annotation" %in% names(model$comp[[1]])) {
    mapping <- try(getAnnMap("SYMBOL", model$comp[[1]]$annotation), silent=TRUE)
    if (class(mapping) != "try-error")
      nameMapping <- mapping
  }

  is_gpdisim_model <- (model$type == 'cgpdisim')
  
  if (is_gpdisim_model) {
    numGenes <- model$comp[[1]]$numGenes + 1
  }
  else {
    numGenes <- model$comp[[1]]$numGenes
  }

  genes <- model$comp[[1]]$genes
  tf <- genes[1]
  targetGenes <- genes[2:length(genes)]

  if (numGenes <= 2) {
    v <- matrix(seq(1, (numGenes+1)*length(model$comp)), nrow=(numGenes+1))
    v[1:2,] <- v[seq(2, 1, by=-1),]
    layout(v)
    singlePlot <- TRUE
  }
  else {
    singlePlot <- FALSE
  }
  model <- modelUpdateProcesses(model, predt=predt)
  for ( i in seq(along=model$comp) ) {
    while (any(model$comp[[i]]$varF < 0) || any(model$comp[[i]]$ypredVar < 0)) {
      warning('Negative variances in GPPlot, adding jitter to fix')
      if (any(model$comp[[i]]$varF < 0))
        model$comp[[i]]$varF <- model$comp[[i]]$varF - 2*min(model$comp[[i]]$varF)
      if (any(model$comp[[i]]$ypredVar < 0))
        model$comp[[i]]$ypredVar <- model$comp[[i]]$ypredVar - 2*min(model$comp[[i]]$ypredVar)
    }

    if (!singlePlot)
      par(mfrow = c(3, ceiling((numGenes+1) / 3)))

    #par(mfrow = c(2, trunc(numPlots / 2 + 0.5)))
    #if (!fileOutput)
    #  dev.new()

    #par(mfrow=c(2, 2))
    plot(model$comp[[i]]$predt, model$comp[[i]]$predF,
         ylim=c(min(model$comp[[i]]$predF-2*sqrt(model$comp[[i]]$varF)),
           max(model$comp[[i]]$predF+2*sqrt(model$comp[[i]]$varF))),
         type="l", lwd=3, xlab="Time",ylab="")
    title("Inferred TF Protein Concentration")
    lines(model$comp[[i]]$predt, model$comp[[i]]$predF+2*sqrt(model$comp[[i]]$varF), lty=2, lwd=3, col=2)
    lines(model$comp[[i]]$predt, model$comp[[i]]$predF-2*sqrt(model$comp[[i]]$varF), lty=2, lwd=3, col=2)

    for ( j in seq(length=numGenes) ) {
      plot(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j],
           ylim=c(min(c(model$comp[[i]]$ypred[,j]-2*sqrt(model$comp[[i]]$ypredVar[,j]),
                        model$comp[[i]]$y[,j]-2*sqrt(model$comp[[i]]$yvar[,j]))),
             max(c(model$comp[[i]]$ypred[,j]+2*sqrt(model$comp[[i]]$ypredVar[,j]),
                   model$comp[[i]]$y[,j]+2*sqrt(model$comp[[i]]$yvar[,j])))),
           type="l", lwd=3, xlab="Time",ylab="")
      if ( is.null(nameMapping) ) {
        genename <- genes[j]
      } else {
        genename <- paste(get(genes[j], nameMapping), " (", genes[j], ")", sep='')
      }
      
      if ( j==1 && is_gpdisim_model ) {
        title(paste(genename, "mRNA (input)"))
      } else {
        title(paste(genename, "mRNA"))
      }

      if (is.null(plotTime)) {
        if ("realt" %in% names(model$comp[[i]]))
          plotTime <- model$comp[[i]]$realt
        else if ("timesCell" %in% names(model$comp[[i]]))
          plotTime <- model$comp[[i]]$timesCell$mRNA
        else
          plotTime <- model$comp[[i]]$t
      }

      # plotCI seems to generate a lot of spurious warnings
      warnOption <- getOption('warn')
      options(warn=-1)
      plotCI(plotTime, model$comp[[i]]$y[,j],
             uiw=2*sqrt(model$comp[[i]]$yvar[,j]), lwd=3, col=3, add=TRUE)
      options(warn=warnOption)
      #points(model$comp[[i]]$t, model$comp[[i]]$y[,j], lwd=3, col=3)
      lines(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j]+2*sqrt(model$comp[[i]]$ypredVar[,j]), lty=2, lwd=3, col=2)
      lines(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j]-2*sqrt(model$comp[[i]]$ypredVar[,j]), lty=2, lwd=3, col=2)
    }
  }
 
}


plotTimeseries <- function(data, nameMapping=NULL) {
  require(annotate)
  
  mapping <- try(getAnnMap("SYMBOL", annotation(data)), silent=TRUE)
  if (class(mapping) != "try-error")
    nameMapping <- mapping

  exps <- unique(data$experiments)
  tshifts <- seq(-0.2, 0.2, length.out=length(exps))

  ngenes <- dim(data)[1]

  nrows = floor(sqrt(ngenes))
  par(mfrow = c(nrows, ceiling(ngenes / nrows)))

  for (j in seq(ngenes)) {
    for (i in exps) {
      t <- data$modeltime[data$experiments==i] + tshifts[i]
      m <- exprs(data)[j,data$experiments==i]
      if (! is.null(var.exprs(data))) {
        sd <- sqrt(var.exprs(data)[j,data$experiments==i])
        maxsd <- sqrt(var.exprs(data)[j,])
      } else {
        sd <- 0*m
        maxsd <- 0*exprs(data)[j,]
      }
      probe <- featureNames(data)[j]
      if (i==1) {
        plot(t, m, type='l',
             ylim=c(0, max(exprs(data)[j,] + 2*maxsd)))
      } else {
        lines(t, m)
      }
      warnOption <- getOption('warn')
      options(warn=-1)
      plotCI(t, m,
             uiw=2*sd, lwd=1, col=3, add=TRUE)
      options(warn=warnOption)
      if (!is.null(nameMapping)) {
        title(paste(get(probe, nameMapping), " (", probe, ")", sep=''))
      } else {
        title(probe)
      }
    }
  }
}
