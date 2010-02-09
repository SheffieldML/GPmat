GPplot <- function(data, savepath = '', doParams = FALSE,
                   nameMapping = NULL, predt = NULL, fileOutput=FALSE) {
  require(gplots)
  
  FONTSIZE <- 10;
  LINEWIDTH <- 1;
  MARKERSIZE <- 10;

  if (is.list(data))
    data <- data[[1]]
  
  if (is.GPmodel(data))
    model <- data@model
  else {
    if ("model" %in% names(data))
      model <- data$model
    else
      model <- data
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

  for ( i in seq(along=model$comp) ) {
    if (is_gpdisim_model) {
      model$comp[[i]] <- gpdisimUpdateProcesses(model$comp[[i]], predt=predt)
    }
    else {
      model$comp[[i]] <- gpsimUpdateProcesses(model$comp[[i]])
    }

    if (any(model$comp[[i]]$varF < 0) || any(model$comp[[i]]$ypredVar < 0)) {
      warning('Negative variances in GPplot')
      return()
    }
    #par(mfrow = c(2, trunc(numPlots / 2 + 0.5)))
    if (!fileOutput)
      dev.new()

    par(mfrow=c(2, 2))
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
        genename <- nameMapping[genes[j]][[1]]
      }
      
      if ( j==1 && is_gpdisim_model ) {
        title(paste(genename, "mRNA (input)"))
      } else {
        title(paste(genename, "mRNA"))
      }
      plotCI(model$comp[[i]]$realt, model$comp[[i]]$y[,j],
             uiw=2*sqrt(model$comp[[i]]$yvar[,j]), lwd=3, col=3, add=TRUE)
      #points(model$comp[[i]]$t, model$comp[[i]]$y[,j], lwd=3, col=3)
      lines(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j]+2*sqrt(model$comp[[i]]$ypredVar[,j]), lty=2, lwd=3, col=2)
      lines(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j]-2*sqrt(model$comp[[i]]$ypredVar[,j]), lty=2, lwd=3, col=2)
    }
  }
 
}
