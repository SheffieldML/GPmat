GPplot <- function(data, savepath = '', doParams = FALSE, selectGenes = data$genes) {
  require(gplots)
  
  FONTSIZE <- 10;
  LINEWIDTH <- 1;
  MARKERSIZE <- 10;

  is_gpdisim_model <- (model$type == 'cgpdisim')
  
  model <- data$model

  if (is_gpdisim_model) {
    numGenes <- model$comp[[1]]$numGenes + 1
  }
  else {
    numGenes <- model$comp[[1]]$numGenes
  }

  numPlots <- length(selectGenes)
  genes <- data$genes
  tf <- genes[1]
  targetGenes <- genes[2:length(genes)]

  for ( i in seq(along=model$comp) ) {
    if (is_gpdisim_model) {
      model$comp[[i]] <- gpdisimUpdateProcesses(model$comp[[i]])
    }
    else {
      model$comp[[i]] <- gpsimUpdateProcesses(model$comp[[i]])
    }
    #par(mfrow = c(2, trunc(numPlots / 2 + 0.5)))
    par(mfrow=c(2, 2))
    plot(model$comp[[i]]$predt, model$comp[[i]]$predF,
         ylim=c(min(model$comp[[i]]$predF-2*sqrt(model$comp[[i]]$varF)),
           max(model$comp[[i]]$predF+2*sqrt(model$comp[[i]]$varF))),
         type="l", lwd=3, xlab="Time",ylab="")
    title("Inferred TF Protein Concentration")
    lines(model$comp[[i]]$predt, model$comp[[i]]$predF+2*sqrt(model$comp[[i]]$varF), lty=2, lwd=3, col=2)
    lines(model$comp[[i]]$predt, model$comp[[i]]$predF-2*sqrt(model$comp[[i]]$varF), lty=2, lwd=3, col=2)

    for ( j in seq(length=numGenes) ) {
      #x11()
      plot(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j],
           ylim=c(min(c(model$comp[[i]]$ypred[,j]-2*sqrt(model$comp[[i]]$ypredVar[,j]),
                        model$comp[[i]]$y[,j]-2*sqrt(model$comp[[i]]$yvar[,j]))),
             max(c(model$comp[[i]]$ypred[,j]+2*sqrt(model$comp[[i]]$ypredVar[,j]),
                   model$comp[[i]]$y[,j]+2*sqrt(model$comp[[i]]$yvar[,j])))),
           type="l", lwd=3, xlab="Time",ylab="")
      if ( j==1 && is_gpdisim_model ) {
        title(paste(genes[j], "mRNA (input)"))
      } else {
        title(paste(genes[j], "mRNA"))
      }
      plotCI(model$comp[[i]]$t, model$comp[[i]]$y[,j],
             uiw=2*sqrt(model$comp[[i]]$yvar[,j]), lwd=3, col=3, add=TRUE)
      #points(model$comp[[i]]$t, model$comp[[i]]$y[,j], lwd=3, col=3)
      lines(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j]+2*sqrt(model$comp[[i]]$ypredVar[,j]), lty=2, lwd=3, col=2)
      lines(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j]-2*sqrt(model$comp[[i]]$ypredVar[,j]), lty=2, lwd=3, col=2)
    }
  }
 
}
