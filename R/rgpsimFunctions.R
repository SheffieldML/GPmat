gpsimRank <- function(model, times, ry, ryvar, genes) {
  if (model$type=="cgpsim") {
    Nrep <- length(model$comp)
    Nrank <- dim(ry[[1]])[2]
    Ntime <- dim(ry[[1]])[1]
    
    rank <- c()
    fullModel <- list(type="rgpsim")
    for ( i in seq(length=Nrank) ) {
      cat (i, ": Computing the ranking score for ", genes[i])
      rankModel <- list(type="cgpsim")
      for ( j in seq(length=Nrep) ) {
        Ntrain <- model$comp[[j]]$numGenes
        
        ## some error checks
        ## checking times
    
        origModel <- model$comp[[j]]
        fixParams <- gpsimExtractParam(origModel, 2)
        NorigParams <- length(fixParams[[1]])
        if ( origModel$includeNoise ) {
          tfKernNparams <- origModel$kern$comp[[1]]$comp[[1]]$nParams
        } else {
          tfKernNparams <- origModel$kern$comp[[1]]$nParams
        }
        
        options <- list()
        y <- cbind(origModel$y, ry[[j]][,i])
        yvar <- cbind(matrix(origModel$yvar, Ntime, Ntrain), ryvar[[j]][,i])
    
        options$includeNoise <- origModel$includeNoise
        options$proteinPrior <- origModel$proteinPrior
        options$bTransform <- origModel$bTransform
        options$optimiser <- origModel$optimiser
        
        Ngenes <- Ntrain+1
        Ntf <- 1
        
        rankModel$comp[[j]] <- gpsimCreate(Ngenes, Ntf, times, y, yvar, options)
        params <- gpsimExtractParam(rankModel$comp[[j]], 2)

        Nparams <- length(params[[1]])
        if ( origModel$includeNoise ) {
          fixIndex <- 1:(tfKernNparams+2*Ntrain)
          fixIndex <- c(fixIndex, (tfKernNparams+2*(Ntrain+1)+1):(Nparams-Ntrain-2))
          fixIndex <- c(fixIndex, (Nparams-Ntrain):(Nparams-1))
          fixValues <- fixParams$values
        }
        
        params[[1]][fixIndex] <- fixValues
        rankModel$comp[[j]]$fix <- list(index=fixIndex, value=fixValues)
        rankModel$comp[[j]] <- gpsimExpandParam(rankModel$comp[[j]], params[[1]])
      }

      ## Optimise the model
      optOptions <- optimiDefaultOptions()
      optOptions$maxit <- 300
      optOptions$optimiser <- "SCG"
      
      # optOptions <- optimiDefaultOptions()
      # optOptions$maxit <- 1000
      # optOptions$display <- FALSE
      optModel <- try( modelOptimise(rankModel, optOptions) )
      # optModel <- rankModel

      flag <- 1
      while ( flag<3 ) {
        reRun <- 0
        if ( is.list(optModel) ) {
          if ( optModel$lnSchFail ) reRun <- 1
        } else if ( !is.list(optModel) ) {
          reRun <- 1
        }
        
        if ( reRun==1 ) {
          params <- gpsimExtractParam(rankModel$comp[[1]])
          if ( flag==1 ) {
            params[-rankModel$comp[[1]]$fix$index] <- c(1, 0, 0, runif(1))
          } else {
            params[-rankModel$comp[[1]]$fix$index] <- 4*(runif(4)-0.5)
          }
          for ( ri in seq(length=Nrep) ) {
            rankModel$comp[[ri]] <- gpsimExpandParam(rankModel$comp[[ri]], params)
          }
          optModel0 <- optModel
          optModel <- try( modelOptimise(rankModel, optOptions) )
          if ( is.list(optModel) ) {
            if ( optModel$llscore > optModel0$llscore )
              optModel <- optModel0
          }
        }
        flag <- flag+1
      }
        
      if ( is.list(optModel) ) {
        rankModel <- optModel
        rankModel$trainingLl <- model$llscore
        rankModel$flatLl <- flatModelOptimise(rankModel)
        if ( !is.na(rankModel$llscore) & !is.na(rankModel$flatLl) ) {
          rankModel$normlisedLl <- -rankModel$llscore+rankModel$trainingLl+rankModel$flatLl
        } else {
          rankModel$normlisedLl <- NA
        }
        cat ("\n", i, ": ", genes[i], " likelihood score: ", rankModel$normlisedLl, "\n\n")
      } else {
        rankModel$llscore <- NA
        rankModel$flatLl <- NA
        rankModel$rankerr <- optModel
        rankModel$normlisedLl <- NA
      }
      ## rank <- c(rank, rankModel$llscore)
      fullModel$comp[[i]] <- rankModel
      fullModel$comp[[i]]$testNo <- i
      fullModel$comp[[i]]$geneName <- genes[i]  
    }
  }
  return (fullModel)
}



flatModelOptimise <- function (model, modelOpt=1) {
  flatModel <- list(type="flat")
  Nrep <- length(model$comp)
  T <- dim(model$comp[[1]]$y)[1]

  if ( modelOpt==1 ) {
    for ( i in seq(length=Nrep) ) {
      testNo <- 6
      flatModel$comp[[i]] <- list()
      flatModel$comp[[i]]$y <- as.matrix(model$comp[[i]]$y[, testNo])
      flatModel$comp[[i]]$yvar <- as.matrix(model$comp[[i]]$yvar[, testNo])
      if ( model$comp[[i]]$includeNoise ) {
        flatModel$comp[[i]]$includeNoise <- model$comp[[i]]$includeNoise
      }

      flatModel$comp[[i]]$mu <- 1 # model$comp[[i]]$mu[testNo]
      flatModel$comp[[i]]$numGenes <- 1
      
      flatModel$comp[[i]]$m <- as.vector(model$comp[[i]]$y - flatModel$comp[[i]]$mu)
      
      if ( model$comp[[i]]$includeNoise ) {
        noise <- whiteKernExtractParam(model$comp[[i]]$kern$comp[[2]]$comp[[testNo+1]])
        
        flatModel$comp[[i]]$noiseVec <- noise
        flatModel$comp[[i]]$K <- diag(as.vector(flatModel$comp[[i]]$yvar+noise))
      } else {
        flatModel$comp[[i]]$K <- diag(as.vector(flatModel$comp[[i]]$yvar[, testNo]))
      }
      
      invK <- jitCholInv(flatModel$comp[[i]]$K)
      flatModel$comp[[i]]$invK <- invK$invM
      flatModel$comp[[i]]$logDetK <- 2* sum( log ( diag(invK$chol) ) )
    }
  } else {
    for ( i in seq(length=Nrep) ) {
      flatModel$comp[[i]] <- list()
      flatModel$comp[[i]]$y <- model$comp[[i]]$y
      flatModel$comp[[i]]$yvar <- model$comp[[i]]$yvar
      if ( model$comp[[i]]$includeNoise ) {
        flatModel$comp[[i]]$includeNoise <- model$comp[[i]]$includeNoise
      }
      
      flatModel$comp[[i]]$mu <- model$comp[[i]]$mu
      flatModel$comp[[i]]$numGenes <- model$comp[[i]]$numGenes
      
      muMat <- t( as.matrix(flatModel$comp[[i]]$mu) %*% matrix(1,1,T) )
      flatModel$comp[[i]]$m <- as.vector(model$comp[[i]]$y - muMat)
      
      if ( model$comp[[i]]$includeNoise ) {
        noise <- c()
        for ( j in seq(length=flatModel$comp[[i]]$numGenes) ) {
          noise <- c(noise, whiteKernExtractParam(model$comp[[i]]$kern$comp[[2]]$comp[[j+1]]))
        }
        flatModel$comp[[i]]$noiseVec <- noise
        
        noiseMat <- t( as.matrix(noise) %*% matrix(1,1,T) )
        flatModel$comp[[i]]$K <- diag(as.vector(flatModel$comp[[i]]$yvar + noiseMat))
      } else {
        flatModel$comp[[i]]$K <- diag(as.vector(flatModel$comp[[i]]$yvar))
      }
      
      invK <- jitCholInv(flatModel$comp[[i]]$K)
      flatModel$comp[[i]]$invK <- invK$invM
      flatModel$comp[[i]]$logDetK <- 2* sum( log ( diag(invK$chol) ) )
    }
  }
  
  ## Optimise the model
  optOptions <- optimiDefaultOptions()
  optOptions$maxit <- 1000
  optOptions$display <- FALSE

  optModel <- try( modelOptimise(flatModel, optOptions) )
  flag <- 1
  # to run with different initialisations, set flag<3
  while ( flag<3 ) {
    reRun <- 0
    if ( is.list(optModel) ) {
      if ( optModel$lnSchFail ) reRun <- 1
    } else if ( !is.list(optModel) ) {
      reRun <- 1
    }

    if ( reRun==1 ) {
      for ( i in seq(length=Nrep) ) {
        flatModel$comp[[i]]$mu <- 2*runif(1)
        flatModel$comp[[i]]$m <- as.vector(flatModel$comp[[i]]$y - flatModel$comp[[i]]$mu)
      }
      optModel <- try( modelOptimise(flatModel, optOptions) )
    }
    flag <- flag+1
  }
  
  if ( is.list(optModel) ) {
    ll <- optModel$llscore
  } else if ( modelOpt==1 ) {
    ll <- NA
  }
  
  return (ll)

}


flatObjective <- function (params, model) {
  T <- dim(model$comp[[1]]$y)[1]
  Nrep <- length(model$comp)

  ll <- 0
  for ( i in seq(length=Nrep) ) {
    model$comp[[i]]$mu <- params
    muMat <- t( as.matrix(model$comp[[i]]$mu) %*% matrix(1,1,T) )
    model$comp[[i]]$m <- as.vector(model$comp[[i]]$y - muMat)
    ll <- ll + gpsimLogLikelihood(model$comp[[i]])
  }
  
  return (-ll)
}



flatGradient <- function (params, model) {
  T <- dim(model$comp[[1]]$y)[1]
  Ngenes <- dim(model$comp[[1]]$y)[2]
  Nrep <- length(model$comp)

  grad <- 0
  for ( i in seq(length=Nrep) ) {
    model$comp[[i]]$mu <- params
    muMat <- t( as.matrix(model$comp[[i]]$mu) %*% matrix(1,1,T) )
    model$comp[[i]]$m <- as.vector(model$comp[[i]]$y - muMat)

    rawGrad <- model$comp[[i]]$invK %*% model$comp[[i]]$m
    gradMat <- matrix(rawGrad, T, Ngenes)
    grad <- grad + apply(gradMat, 2, mean)
  }

  grad <- -grad
  return (grad)
}



flatExtractParam <- function (model) {
  Nrep <- length(model$comp)

  params <- 0
  for ( i in seq(length=Nrep) ) {
    params <- params + model$comp[[i]]$mu
  }
  params <- params/Nrep

  return (params)
}



flatExpandParam <- function (model, params) {
  Nrep <- length(model$comp)

  for ( i in seq(length=Nrep) ) {
    model$comp[[i]]$mu <- params
  }

  return (model)
}



flatOptimise <- function (model, options, ...) {
  newOptions <- options
  options <- list()
  for ( i in 1:4 ) 
    options[[i]] <- newOptions[[i]]
  
  if ( any(grep("display",names(newOptions))) ) {
    options$display <- newOptions$display
  }
    
  params <- flatExtractParam(model)

  if ( newOptions$optimiser=="CG" )
    newParams <- CGoptim(params, fn=flatObjective, grad=flatGradient, options, model)

  model <- flatExpandParam(model, newParams$xmin)
  model$llscore <- newParams$objective
  model$lnSchFail <- newParams$lnSchFail
  
  return (model)
}



rgpsimRankSingleProfile <- function (rankModel, genes, expType, expNo, option=1) {

  # locate the positions in the rank model
  geneFullInfo <- c()
  Ntest <- length(rankModel$comp)
  for ( i in seq(length=Ntest) ) 
    geneFullInfo <- c(geneFullInfo, rankModel$comp[[i]]$geneName)
  
  matchInd <- c()
  for ( i in seq(along=genes) ) {
    ind <- grep(genes[i], geneFullInfo)
    index <- c()
    for ( j in seq(length=length(ind)) ) {
      if ( genes[i] %in% geneFullInfo[ind[j]] )
        index <- c(index, ind[j])      
    }
    if ( length(index) != 1 ) {
      warnmsg <- paste("Too many or too few matches for ",genes[i], "!", sep="")
      warning(warnmsg)
    } else {
      matchInd <- c(matchInd, index)
    }
  }

  if ( option!=1 ) 
    postscript(paste(expType, expNo, ".ps", sep=""), horizontal=FALSE, width=8.0, height=6.0)
  
  for ( i in seq(along=matchInd) ) {
    model <- ranks$comp[[matchInd[i]]]
    for ( i in seq(length=Nrep) ) 
      model$comp[[i]] <- gpsimUpdateProcesses(model$comp[[i]])
    
    ## Display the result
    if ( option==1 ) {
      for ( i in seq(along=model$comp) ) {
        x11()
        ymax <- max(model$comp[[i]]$predF+2*sqrt(model$comp[[i]]$varF))
        ymin <- min(model$comp[[i]]$predF-2*sqrt(model$comp[[i]]$varF))
        plot(model$comp[[i]]$predt, model$comp[[i]]$predF, ylim=c(ymin, ymax), type="l", lwd=3, xlab="Time",ylab="")
        title(paste("Predicted Protein Concentration (for", model$geneName, ")", sep=" "))
        lines(model$comp[[i]]$predt, model$comp[[i]]$predF+2*sqrt(model$comp[[i]]$varF), lty=2, lwd=3, col=2)
        lines(model$comp[[i]]$predt, model$comp[[i]]$predF-2*sqrt(model$comp[[i]]$varF), lty=2, lwd=3, col=2)
        
        j <- 6
        x11()
        ymax <- max(model$comp[[i]]$ypred[,j]+2*sqrt(model$comp[[i]]$ypredVar[,j]))
        ymin <- min(model$comp[[i]]$ypred[,j]-2*sqrt(model$comp[[i]]$ypredVar[,j]))        
        plot(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j], ylim=c(ymin, ymax), type="l", lwd=3, xlab="Time",ylab="")
        title(paste("mRNA", model$geneName))
        points(model$comp[[i]]$timesCell[[j+1]], model$comp[[i]]$y[,j], lwd=3, col=3) 
        lines(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j]+2*sqrt(model$comp[[i]]$ypredVar[,j]), lty=2, lwd=3, col=2)
        lines(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j]-2*sqrt(model$comp[[i]]$ypredVar[,j]), lty=2, lwd=3, col=2)
      }
      
    } else {
      for ( i in seq(along=model$comp) ) {
        ymax <- max(model$comp[[i]]$predF+2*sqrt(model$comp[[i]]$varF))
        ymin <- min(model$comp[[i]]$predF-2*sqrt(model$comp[[i]]$varF))
        plot(model$comp[[i]]$predt, model$comp[[i]]$predF, ylim=c(ymin, ymax), type="l", lwd=3, xlab="Time",ylab="")
        title(paste("Predicted Protein Concentration (for", model$geneName, ")", sep=" "))
        lines(model$comp[[i]]$predt, model$comp[[i]]$predF+2*sqrt(model$comp[[i]]$varF), lty=2, lwd=3, col=2)
        lines(model$comp[[i]]$predt, model$comp[[i]]$predF-2*sqrt(model$comp[[i]]$varF), lty=2, lwd=3, col=2)
        
        j <- 6
        ymax <- max(model$comp[[i]]$ypred[,j]+2*sqrt(model$comp[[i]]$ypredVar[,j]))
        ymin <- min(model$comp[[i]]$ypred[,j]-2*sqrt(model$comp[[i]]$ypredVar[,j]))        
        plot(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j], ylim=c(ymin, ymax), type="l", lwd=3, xlab="Time",ylab="")
        title(paste("mRNA", model$geneName))
        points(model$comp[[i]]$timesCell[[j+1]], model$comp[[i]]$y[,j], lwd=3, col=3) 
        lines(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j]+2*sqrt(model$comp[[i]]$ypredVar[,j]), lty=2, lwd=3, col=2)
        lines(model$comp[[i]]$predt, model$comp[[i]]$ypred[,j]-2*sqrt(model$comp[[i]]$ypredVar[,j]), lty=2, lwd=3, col=2)
      }
    }
  }
  
  if ( option!=1) dev.off()
}



 

