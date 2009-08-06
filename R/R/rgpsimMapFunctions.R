rgpsimMapRank <- function(model, times, ry, ryvar, genes) {
  if (model$type=="cgpsimMap") {
    Nrep <- length(model$comp)
    Nrank <- dim(ry[[1]])[2]
    Ntime <- dim(ry[[1]])[1]
    
    rank <- c()
    fullModel <- list(type="rgpsimMap")
    for ( i in seq(length=Nrank) ) {
      cat (i, ": Computing the ranking score for ", genes[i], "\n")
      rankModel <- list(type="cgpsimMap")
      for ( j in seq(length=Nrep) ) {
        Ntrain <- model$comp[[j]]$numGenes
        
        ## some error checks
        ## checking times
    
        origModel <- model$comp[[j]]
        fixParams <- gpsimMapExtractParam(origModel, 2)
        Nfixparams <- origModel$numParams
        
        options <- list()
        y <- cbind(origModel$y, ry[[j]][,i])
        yvar <- cbind(matrix(origModel$yvar, Ntime, Ntrain), ryvar[[j]][,i])
    
        options$includeNoise <- origModel$includeNoise
        options$proteinPrior <- origModel$proteinPrior
        options$Transform <- origModel$Transform
        options$optimiser <- origModel$optimiser
        options$nonLinearity <- origModel$nonLinearity
        options$step <- origModel$step
        options$mapt <- origModel$mapt
        options$times <- times
        options$timesIndex <- origModel$timesIndex
        
        Ngenes <- Ntrain+1
        Ntf <- 1

        options$S <- array(1, Ngenes)
        set.seed(1)
        options$D <- runif(Ngenes)
        mu <- apply(y, 2, mean)
        options$B <- options$D*mu
        
        options$kern <- origModel$kern
        
        rankModel$comp[[j]] <- gpsimMapCreate(Ngenes, Ntf, times, y, yvar, options)
        rankModel$comp[[j]]$updateF <- FALSE

        params <- gpsimMapExtractParam(rankModel$comp[[j]], 2)
        
        Nparams <- length(params[[1]])
        if ( origModel$includeNoise ) {
          fixIndex <- 1:Nfixparams
          fixValues <- fixParams$values[1:Nfixparams]
        }
        
        params[[1]][fixIndex] <- fixValues
        rankModel$comp[[j]]$fix <- list(index=fixIndex, value=fixValues)
        rankModel$comp[[j]] <- gpsimMapExpandParam(rankModel$comp[[j]], params[[1]])
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
          if ( any(names(optModel)=="lnSchFail") )
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
        rankModel$flatLl <- flatMapOptimise(rankModel)
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
    return (fullModel)
  } else {
    stop("Wrong model type!")
  }

}



flatMapOptimise <- function (model, modelOpt=1) {
  flatModel <- list(type="flatMap")
  Nrep <- length(model$comp)
  T <- dim(model$comp[[1]]$y)[1]

  if ( modelOpt==1 ) {
    for ( i in seq(length=Nrep) ) {
      testNo <- model$comp[[1]]$numGenes
      flatModel$comp[[i]] <- list()
      flatModel$comp[[i]]$y <- as.matrix(model$comp[[i]]$y[, testNo])
      flatModel$comp[[i]]$yvar <- as.matrix(model$comp[[i]]$yvar[, testNo])
      if ( model$comp[[i]]$includeNoise ) {
        flatModel$comp[[i]]$includeNoise <- model$comp[[i]]$includeNoise
      }

      flatModel$comp[[i]]$mu <- 1 # model$comp[[i]]$mu[testNo]
      flatModel$comp[[i]]$numGenes <- 1
      
      flatModel$comp[[i]]$m <- as.vector(model$comp[[i]]$y - flatModel$comp[[i]]$mu)
      browser()

      if ( model$comp[[i]]$includeNoise ) {
        noise <- model$comp[[i]]$noiseVar[testNo]
        flatModel$comp[[i]]$K <- diag(as.vector(flatModel$comp[[i]]$yvar+noise))
      } else {
        flatModel$comp[[i]]$K <- diag(as.vector(flatModel$comp[[i]]$yvar[, testNo]))
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
