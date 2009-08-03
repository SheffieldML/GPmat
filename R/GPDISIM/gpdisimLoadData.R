gpdisimLoadMef2Data <- function () {
  ## locate data folder in the work directory.
  dataPath <- paste(getwd(), "data", sep="/")
  
  fileName <- paste(dataPath, "mef2Data.RData", sep="/")
  if (file.exists(fileName)) {
    load(fileName)
  } else {
    fileName <- paste(dataPath, "mef2Data_y.txt", sep="/")
    yFull <- read.table(file=fileName)

    fileName <- paste(dataPath, "mef2Data_yvar.txt", sep="/")
    yFullVar <- read.table(file=fileName)

    fileName <- paste(dataPath, "mef2Data_scale.txt", sep="/")
    scale <- read.table(file=fileName)

    times <- seq(1, 12)

    # gene names need to be inserted!
    genes <- c()

    y <- list()
    yvar <- list()
    numRep <- 3
    numData <- length(times)
    for ( i in seq(length=numRep) ) {
      y[[i]] <- yFull[((i-1)*numData+1):(i*numData),]
      yvar[[i]] <- yFullVar[((i-1)*numData+1):(i*numData),]
    }

    trainingData <- list(y=y, yvar=yvar, genes=genes, times=times, scale=scale)

    pwd <- getwd()
    setwd(dataPath)
    save(trainingData, file="mef2Data.RData")
    setwd(pwd)
  }

  return (trainingData)
}
