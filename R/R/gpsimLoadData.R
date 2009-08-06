gpsimLoadEcoliFullData <- function () {
  ## locate data folder in the work directory.
  dataPath <- paste(getwd(), "data", sep="/")
  
  fileName <- paste(dataPath, "ecoliNormalisedData.txt", sep="/")
  if (file.exists(fileName)) {
    ecoliData <- read.table(file=fileName)
    geneNames <- c("dinF", "dinI", "lexA", "recA", "recN", "ruvA", "ruvB", "sbmC", "sulA", "umuC", "umuD", "uvrB", "yebG", "yjiW")
    rawExp <- list(data=t(ecoliData), genes=geneNames)
    
    targetInd <- 1:14
    yFull <- rawExp$data[,targetInd]

    Ngenes <- max(targetInd)
    Nrep <- 1
    Nt <- length(yFull[,1])

    ## Rescale so that the SD of the curves is 1
    scale <- sqrt(apply(yFull, 2, var))
    scaleMat <- t( matrix(scale, Ngenes, Nrep*Nt) )
    y <- list()
    y[[1]] <- yFull/scaleMat
    times <- c(0, 5, 10, 20, 40, 60)
    ## times <- times/2
    genes <- rawExp$genes[targetInd]

  } else {	
    stop("The source data file cannot be found!")
  }

  trainingData <- list(y=y, genes=genes, times=times, scale=scale)
  return (trainingData)
}





gpsimLoadBarencoPUMAData <- function(option=1) {
  ## locate data folder in the work directory.
  dataPath <- paste(getwd(), "data","p53BarencoData", sep="/")

  if (option==1) {
    fileName <- paste(dataPath, "barencoPUMAtraining.RData", sep="/")
  } else {
    fileName <- paste(dataPath, "barencoPUMAfull.RData", sep="/")
  }

  if (file.exists(fileName)) {
    load(fileName)
    if (option==1) {
      return (trainingData)
    } else {
      return (importData)
    } 
  } else {
    if (option==1) {
      fileName <- paste(dataPath, "p53BarencoPUMA.txt", sep="/")
      if ( file.exists(fileName) ) {
        oldPath <- getwd()
        setwd(dataPath)
        p53Data <- read.table(file="p53BarencoPUMA.txt")
        trainingGenes <- c("203409_at", "205780_at", "209295_at", "202284_s_at", "218346_s_at")
        numGenes <- 5
        numData <- 7
        numRep <- 3      
        fullExprs <- matrix(0, numData*numRep, numGenes)
        fullVar <- matrix(0, numData*numRep, numGenes)
        scale <- array(1, numGenes)
        for ( k in seq(length=numGenes) ) {
          fullExprs[,k] <- p53Data[[k+1]][1:(numData*numRep)]
          fullVar[,k] <- p53Data[[k+1]][(numData*numRep+1):(numData*numRep*2)]
          scale[k] <- 1e3*p53Data[43, k+1]
        }
        
        y <- list()
        yvar <- list()
        for ( i in seq(length=numRep) ) {
          y[[i]] <- fullExprs[((i-1)*numData+1):(i*numData),]
          yvar[[i]] <- fullVar[((i-1)*numData+1):(i*numData),]
        }
        
        genes <- list()      
        genes$tag <- c("203409_at", "205780_at", "209295_at", "202284_s_at", "218346_s_at")
        genes$names <- c("DDB2", "BIK", "TNFRSF10b", "p21", "hPA26")
        
        times <- seq(0, 12, by=2)
        
        trainingData <- list(y=y, yvar=yvar, genes=genes, times=times, scale=scale)
            
        save(trainingData, file="barencoPUMAtraining.RData")
        setwd(oldPath)
      
        return (trainingData)
      } else {
        stop("No source file is found!")
      }
    } else if (option==2) {
      csvFiles <- list.files(path=dataPath, pattern="csv")
      fileNames <- c("barencoPUMA_exprs.csv", "barencoPUMA_se.csv")
      ind <- match(fileNames, csvFiles)
      if ( any(is.na(ind)) ) {
        library(puma)
        library(Biobase)
        celpath <- paste(dataPath, "barencoRawData", sep="/")
        fns <- list.celfiles(path=celpath, full.names=TRUE)
        cat("loading .CEL files ... \n")
        abatch <- ReadAffy(celfile.path=celpath)
        
        eset_mmgmos <- mmgmos(abatch)
        write.reslts(eset_mmgmos, file="barencoPUMA")
      } 
      
      fileName <- paste(dataPath, "barencoPUMA_exprs.csv", sep="/")
      importExprs <- read.csv(fileName, header=TRUE)
      fileName <- paste(dataPath, "barencoPUMA_se.csv", sep="/")    
      importSd <- read.csv(fileName, header=TRUE)
      
      fileName <- paste(dataPath, "E-MEXP-549-2columns.csv", sep="/")
      columnNames <- read.csv(fileName, header=TRUE)    
      
      columnSortNames <- columnNames[order(columnNames$Index),]
      celSortNames <- columnSortNames$Array.Data.File
      columns <- colnames(importExprs)
      
      fullExprs <- matrix(0, dim(importExprs)[1], dim(importExprs)[2]-1)
      fullSd <- matrix(0, dim(importSd)[1], dim(importSd)[2]-1)
      for ( i in seq(along=celSortNames) ) {
        ind <- match(celSortNames[i], columns)
        fullExprs[,i] <- importExprs[,ind]
        fullSd[,i] <- importSd[,ind]      
      }
      
      fullExprs <- t(fullExprs)
      fullSd <- t(fullSd)
      fullVar <- fullSd*fullSd
      
      genes <- importExprs$X
      
      yFull <- exp(fullExprs + fullVar/2)
      yVarFull <- (exp(fullVar)-1)*exp(2*fullExprs + fullVar)
      
      ## rescale
      scale <- sqrt( apply(yFull,2,var) )
      Ngene <- dim(yFull)[2]
      Nrep <- 3
      Nt <- 7
      scaleMat <- t( matrix(scale, Ngene, Nrep*Nt) )
      yFull <- yFull/scaleMat
      yVarFull <- yVarFull/(scaleMat*scaleMat)

      y <- list()
      yVar <- list()
      for ( i in 1:Nrep ) {
        y[[i]] <- yFull[(i-1)*Nt+1:Nt,]
        yVar[[i]] <- yVarFull[(i-1)*Nt+1:Nt,]
      }

      times <- seq(0, 12, by=2)
      
      importData <- list(y=y, yvar=yVar, genes=genes, times=times, scale=scale)
      
      pwd <- getwd()
      setwd(dataPath)
      save(importData, file="barencoPUMAfull.RData")
      setwd(pwd)

      return (importData)      
    }
  }
}


gpsimLoadBarencoMAS5Data <- function(option=1) {
  ## option = 1: load training data
  ## option = 2: load target ranking data
  ## option = 3: load full data set
  ## locate data folder in the work directory.
  dataPath <- paste(getwd(), "data","p53BarencoData", sep="/")

  if (option==1) {
    fileName <- paste(dataPath, "barencoMAS5training.RData", sep="/")
  } else if (option==2) {
    fileName <- paste(dataPath, "barencoMAS5ranking.RData", sep="/")
  } else {
    fileName <- paste(dataPath, "barencoMAS5full.RData", sep="/")
  }
  
  if (file.exists(fileName)) {
    load(fileName)
    return (importData)
  } else {
      fileName <- paste(dataPath, "barencoMAS5.RData", sep="/")
      if (!file.exists(fileName))
        stop("The source data file cannot be found!") else
      ## load experiment source data which is an ExpressionSet object 
      load(fileName)

      fullGeneNames <- featureNames(fiveGyMAS5)
      NtrainGenes <- 5
      matchInd <- numeric(0)
      ## Select the training genes
      for ( i in 1:NtrainGenes ) {
        ind <- grep(genes$tag[i], fullGeneNames)
        index <- c()
        for ( j in seq(length=length(ind)) ) {
          if ( genes$tag[i] %in% fullGeneNames[ind[j]] )
            index <- c(index, ind[j])      
        }
        if ( length(index) != 1 )
          stop("Too many or too few matches!") else
        matchInd <- c(matchInd, index)
      }

      fullExprs <- assayDataElement(fiveGyMAS5, "exprs")
      fullSd <- assayDataElement(fiveGyMAS5, "se.exprs")
      
      phenoData <- pData(fiveGyMAS5)
      Nrep <- max(phenoData$replicate)
      
      trainExprs <- fullExprs[matchInd,]
      trainSd <- fullSd[matchInd,]

      if (option==1) {
        exprsData <- trainExprs
        sdData <- trainSd
        ## number of the training genes
        Ngene <- length(genes$names)
      } else if (option==2) {
        exprsData <- fullExprs[-matchInd,]
        sdData <- fullSd[-matchInd,]
        genes <- attr(exprsData, "dimnames")[[1]]
        Ngene <- length(genes)
      } else {
        exprsData <- fullExprs
        sdData <- fullSd
        genes <- attr(exprsData, "dimnames")[[1]]
        Ngene <- length(genes)
      }
      
      y <- list()
      yVar <- list()
      for ( i in 1:Nrep ) {
        repInd <- grep(i, phenoData$replicate)

        celFiles <- attr(phenoData, "row.names")[repInd]

        timeInfo <- sort(phenoData$time[repInd], index.return=TRUE)

        times <- timeInfo$x
        Nt <- length(times)
        yRepi <- list(matrix(0, Nt, Ngene))
        yVarRepi <- list(matrix(0, Nt, Ngene))

        for ( ti in 1:Nt ) {
          exprsCel <- attr(exprsData, "dimnames")[[2]]
          exprsInd <- grep(celFiles[timeInfo$ix[ti]], exprsCel)

          yRepi[[1]][ti,] <- exprsData[,exprsInd]
          yVarRepi[[1]][ti,] <- sdData[,exprsInd]
        }

        y <- c(y, yRepi)
        yVar <- c(yVar, yVarRepi)
      }

      yFull <- y[[1]]
      yVarFull <- yVar[[1]]
      for ( i in 2:Nrep ) {
        yFull <- rbind(yFull, y[[i]])
        yVarFull <- rbind(yVarFull, y[[i]])
      }

      ## rescale
      scale <- sqrt( apply(yFull,2,var) )
      scaleMat <- t( matrix(scale, Ngene, Nrep*Nt) )
      yFull <- yFull/scaleMat
      yVarFull <- yVarFull/(scaleMat*scaleMat)

      for ( i in 1:Nrep ) {
        y[[i]] <- yFull[(i-1)*Nt+1:Nt,]
        yVar[[i]] <- yVarFull[(i-1)*Nt+1:Nt,]
      }      

      importData <- list(y=y, yvar=yVar, genes=genes, times=times, scale=scale)

      pwd <- getwd()
      setwd(dataPath)
      
      if (option==1) {
        save(importData, file="barencoMAS5training.RData")
      } else if (option==2) {
        save(importData, file="barencoMAS5ranking.RData")
      } else {
        save(importData, file="barencoMAS5full.RData")
      }      

      setwd(pwd)
      
      return (importData)
    }
}
  


gpsimLoadElkData <- function(option=1) {
  ## locate data folder in the work directory.
  dataPath <- paste(getwd(), "data","ResultsAmitHela_mmgmos", sep="/")

  if (option==1) {
    fileName <- paste(dataPath, "elktraining1.RData", sep="/")
  } else  if (option==2) {
    fileName <- paste(dataPath, "elktraining2.RData", sep="/")
  } else  if (option==3) {
    fileName <- paste(dataPath, "elktraining3.RData", sep="/")
  } else {
    fileName <- paste(dataPath, "elkfull.RData", sep="/")
  }

  if (file.exists(fileName)) {
    load(fileName)
    if (option<=3) {
      return (trainingData)
    } else {
      return (importData)
    } 
  } else {
    csvFiles <- list.files(path=dataPath, pattern="csv")
    fileNames <- c("resultsAmitHela_exprs.csv", "resultsAmitHela_se.csv")
    ind <- match(fileNames, csvFiles)
    if ( any(is.na(ind)) ) 
      stop("No source csv files are found!")
    
    fileName <- paste(dataPath, "resultsAmitHela_exprs.csv", sep="/")
    importExprs <- read.csv(fileName, header=TRUE)
    fileName <- paste(dataPath, "resultsAmitHela_se.csv", sep="/")    
    importSd <- read.csv(fileName, header=TRUE)

    Ngenes <- dim(importExprs)[1]
    Nt <- dim(importExprs)[2]-1

    fullExprs <- importExprs[,2:(Nt+1)]
    fullSd <- importSd[,2:(Nt+1)]
      
    fullExprs <- t(fullExprs)
    fullSd <- t(fullSd)
    fullVar <- fullSd*fullSd
      
    genes <- importExprs$X
    
    yFull <- exp(fullExprs + fullVar/2)
    yVarFull <- (exp(fullVar)-1)*exp(2*fullExprs + fullVar)
    
    ## rescale
    scale <- sqrt( apply(yFull,2,var) )
    Ngene <- dim(yFull)[2]
    Nrep <- 1
    scaleMat <- t( matrix(scale, Ngene, Nrep*Nt) )
    yFull <- yFull/scaleMat
    yVarFull <- yVarFull/(scaleMat*scaleMat)

    times <- c(0, 1/3, 2/3, 1, 2, 4, 8)

    if ( option<=3 ) {
      if ( option == 1 ) {
        trainingGenes1 <- c("209101_at", "201280_s_at", "204897_at", "201617_x_at", "210764_s_at", "208961_s_at")
      } else if ( option == 2 ) {
        trainingGenes1 <- c("212721_at", "219622_at", "202431_s_at", "211090_s_at", "219793_at", "219130_at")
      } else {
        trainingGenes1 <- c("201693_s_at", "209189_at", "204420_at", "202081_at", "201473_at", "200797_s_at", "202340_x_at")
      }
      Ntrain <- length(trainingGenes1)
      matchInd <- c()
      for ( i in 1:Ntrain ) {
        ind <- grep(trainingGenes1[i], genes)
        index <- c()
        for ( j in seq(length=length(ind)) ) {
          if ( trainingGenes1[i] %in% genes[ind[j]] )
            index <- c(index, ind[j])      
        }
        if ( length(index) != 1 ) {
          warnmsg <- paste("Too many or too few matches for ", affInd[i], "!", sep="")
          warning(warnmsg)
        } else {
          matchInd <- c(matchInd, index)
        }
      }

      y <- list()
      yVar <- list()
      for ( i in 1:Nrep ) {
        y[[i]] <- yFull[(i-1)*Nt+1:Nt, matchInd]
        yVar[[i]] <- yVarFull[(i-1)*Nt+1:Nt, matchInd] 
      }

      scale <- scale[matchInd]
      trainingData <- list(y=y, yvar=yVar, genes=trainingGenes1, times=times, scale=scale)
      pwd <- getwd()
      setwd(dataPath)

      if (option==1) {
        fileName <- paste(dataPath, "elktraining1.RData", sep="/")
      } else  if (option==2) {
        fileName <- paste(dataPath, "elktraining2.RData", sep="/")
      } else  if (option==3) {
        fileName <- paste(dataPath, "elktraining3.RData", sep="/")
      }
      
      save(trainingData, file=fileName)
      setwd(pwd)

      return (trainingData)
      
    } else {
      y <- list()
      yVar <- list()
      for ( i in 1:Nrep ) {
        y[[i]] <- yFull[(i-1)*Nt+1:Nt,]
        yVar[[i]] <- yVarFull[(i-1)*Nt+1:Nt,] 
      }
    
      importData <- list(y=y, yvar=yVar, genes=genes, times=times, scale=scale)
      
      pwd <- getwd()
      setwd(dataPath)
      save(importData, file="elkfull.RData")
      setwd(pwd)
      
      return (importData)      
    }
  }
}
