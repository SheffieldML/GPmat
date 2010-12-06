gpLoadToolkit <- function(gppath=NULL) {
  library(spam)
  #library(tigre)

  if (is.null(gppath)) {
    if (.Platform$OS.type == "unix") {
      dir = "~/mlprojects/gp/R/R/"
    } else if (.Platform$OS.type == "windows") {
      dir = "C:\\mlprojects\\gp\\R\\R\\"
    }
  } else dir = gppath

  ## source() a bunch of files
  sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
      if(trace) cat(nm,":")           
      source(file.path(path, nm), ...)
      if(trace) cat("\n")
    }
  }

  ## Update symbolic links.
#   system("ln -fs ~/mlprojects/gpsim/R/R/cmpndKernFunctions.R ~/mlprojects/gp/R/cmpndKernFunctions.R")
#   system("ln -fs ~/mlprojects/gpsim/R/R/kernFunctions.R ~/mlprojects/gp/R/kernFunctions.R")
#   system("ln -fs ~/mlprojects/gpsim/R/R/miscFunctions.R ~/mlprojects/gp/R/miscFunctions.R")
#   system("ln -fs ~/mlprojects/gpsim/R/R/optimiFunctions.R ~/mlprojects/gp/R/optimiFunctions.R")
#   system("ln -fs ~/mlprojects/gpsim/R/R/rbfKernFunctions.R ~/mlprojects/gp/R/rbfKernFunctions.R")
#   system("ln -fs ~/mlprojects/gpsim/R/R/whiteKernFunctions.R ~/mlprojects/gp/R/whiteKernFunctions.R")

  sourceDir(dir)
}