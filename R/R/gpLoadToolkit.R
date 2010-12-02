gpLoadToolkit <- function() {
  library(spam)
  #library(tigre)

  if (.Platform$OS.type == "unix") {
    dir = "~/mlprojects/gp/R/R/"
  } else if (.Platform$OS.type == "windows") {
    dir = "C:\\mlprojects\\gp\\R\\R\\"
  }

  ## source() a bunch of files
  sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
      if(trace) cat(nm,":")           
      source(file.path(path, nm), ...)
      if(trace) cat("\n")
    }
  }

  ## Update symbolic links.
  #system("ln -fs ~/mlprojects/gpsim/R/R/AllClasses.R ~/mlprojects/gp/R/AllClasses.R")
  system("ln -fs ~/mlprojects/gpsim/R/R/cmpndKernFunctions.R ~/mlprojects/gp/R/cmpndKernFunctions.R")
  system("ln -fs ~/mlprojects/gpsim/R/R/kernFunctions.R ~/mlprojects/gp/R/kernFunctions.R")
  system("ln -fs ~/mlprojects/gpsim/R/R/miscFunctions.R ~/mlprojects/gp/R/miscFunctions.R")
  system("ln -fs ~/mlprojects/gpsim/R/R/optimiFunctions.R ~/mlprojects/gp/R/optimiFunctions.R")
  system("ln -fs ~/mlprojects/gpsim/R/R/rbfKernFunctions.R ~/mlprojects/gp/R/rbfKernFunctions.R")
  system("ln -fs ~/mlprojects/gpsim/R/R/whiteKernFunctions.R ~/mlprojects/gp/R/whiteKernFunctions.R")

  sourceDir(dir)

#   source(paste(dir,"gpOptions.R",sep=""))
#   source(paste(dir,"gpCreate.R",sep=""))
#   source(paste(dir,"gpExtractParam.R",sep=""))
#   source(paste(dir,"gpExpandParam.R",sep=""))
#   source(paste(dir,"miscFunctions.R",sep=""))
#   source(paste(dir,"kernFunctions.R",sep=""))
#   source(paste(dir,"cmpndKernFunctions.R",sep=""))
#   source(paste(dir,"rbfKernFunctions.R",sep=""))
#   source(paste(dir,"whiteKernFunctions.R",sep=""))
#   source(paste(dir,"gpUpdateKernels.R",sep=""))
#   source(paste(dir,"gpUpdateAD.R",sep=""))
#   source(paste(dir,"gpComputeM.R",sep=""))
#   source(paste(dir,"gpOut.R",sep=""))
#   source(paste(dir,"gpPosteriorMeanVar.R",sep=""))
#   source(paste(dir,"gpComputeAlpha.R",sep=""))
#   source(paste(dir,"gpDataIndices.R",sep=""))
#   source(paste(dir,"gpBlockIndices.R",sep=""))
#   source(paste(dir,"gpOptimise.R",sep=""))
#   source(paste(dir,"gpObjective.R",sep=""))
#   source(paste(dir,"gpLogLikelihood.R",sep=""))
#   source(paste(dir,"gpGradient.R",sep=""))
#   source(paste(dir,"gpLogLikeGradients.R",sep=""))
#   source(paste(dir,"gpScaleBiasGradient.R",sep=""))
#   source(paste(dir,"gpTest.R",sep=""))
#   source(paste(dir,"gpCovGradsTest.R",sep=""))
#   source(paste(dir,"gpCovGrads.R",sep=""))
#   source(paste(dir,"modelGradientCheck.R",sep=""))
#   source(paste(dir,"rbfKernGradX.R",sep=""))
#   source(paste(dir,"whiteKernGradX.R",sep=""))
#   source(paste(dir,"rbfKernDiagGradX.R",sep=""))
#   source(paste(dir,"whiteKernDiagGradX.R",sep=""))
#   source(paste(dir,"kernDiagGradient.R",sep=""))
#   source(paste(dir,"gpPlot.R",sep=""))
#   source(paste(dir,"gaussSamp.R",sep=""))
#   source(paste(dir,"zeroAxes.R",sep=""))

}