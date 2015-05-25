#
# selproj kernel
#   wrap around another simple kernel
#   returning its value if x[1,]==kern$options$expmask,
#   or x[1,]==0 or kern$options$expmask==0
#


selprojKernParamInit <- function (kern) {
  kern <- cmpndKernParamInit(kern)

  for (i in seq_along(kern$comp))
    if (kern$comp[[i]]$type == "multi")
      error("selproj kernel does not support wrapping a multi kernel.")

  stopifnot(!is.null(kern$options$expmask))
  
  kern$expmask <- kern$options$expmask
  kern$masklen <- length(kern$expmask)
  kern$transforms <- kern$comp[[1]]$transforms
  if ("transformArgs" %in% names(kern$comp[[1]]))
    kern$transformArgs <- kern$comp[[1]]$transformArgs
  
  return (kern)
}


.selprojDataMaskCompare <- function (m, m1) {
  m1 <- as.matrix(m1)
  I <- vector(len=dim(m1)[1])
  I[] <- TRUE

  for (i in seq_along(m)) {
    if (m[i] > 0)
      I <- I & (m1[,i] == m)
  }
  return (I)
}


.selprojKernMaskCombine <- function (m1, m2) {
  stopifnot(length(m1) == length(m2))

  if (any((m1 > 0) & (m2 > 0) & (m1 != m2)))
    return (NA)
  else
    return (pmax(m1, m2))
}



selprojKernDiagCompute <- function (kern, x) {
  d <- dim(as.array(x))[1]
  k <- array(0, dim=d)

  if (dim(as.array(x))[2] > 1) {
    I <- .selprojDataMaskCompare(kern$expmask, x[,seq(kern$masklen)])

    k[I] <- kernDiagCompute(kern$comp[[1]], x[I,-seq(kern$masklen)])
  }

  return (k)
}



selprojKernExtractParam <- function (kern, only.values=TRUE,
                                     untransformed.values=FALSE) {
  funcName <- paste(kern$comp[[1]]$type, "KernExtractParam", sep="")
  func <- get(funcName, mode="function")

  params <- func(kern$comp[[1]], only.values=only.values,
                 untransformed.values=untransformed.values)

  if (!only.values) {
    names(params) <- paste(kern$comp[[1]]$type, names(params), sep='_')
  }

  return (params)
}



selprojKernExpandParam <- function (kern, params) {
  funcName <- paste(kern$comp[[1]]$type, "KernExpandParam", sep="")
  func <- get(funcName, mode="function")
  kern$comp[[1]] <- func(kern$comp[[1]], params)
  
  return (kern)
}



selprojKernCompute <- function (kern, x, x2) {
  if ( nargs()>2 ) {
    d1 <- dim(as.array(x))[1]
    d2 <- dim(as.array(x2))[1]
    k <- array(0, dim=c(d1, d2))

    if (dim(as.array(x))[2] > 1 && dim(as.array(x2))[2] > 1) {
      I1 <- .selprojDataMaskCompare(kern$expmask, x[,seq(kern$masklen)])
      I2 <- .selprojDataMaskCompare(kern$expmask, x2[,seq(kern$masklen)])

      if (!any(I1) || !any(I2))
        return (k)

      k[I1,I2] <- kernCompute(kern$comp[[1]],
                              x[I1,-seq(kern$masklen)],
                              x2[I2,-seq(kern$masklen)])
    }
  } else {
    d1 <- dim(as.array(x))[1]
    k <- array(0, dim=c(d1, d1))

    if (dim(as.array(x))[2] > 1) {
      I1 <- .selprojDataMaskCompare(kern$expmask, x[,seq(kern$masklen)])

      if (!any(I1))
        return (k)
      
      k[I1,I1] <- kernCompute(kern$comp[[1]], x[I1,-seq(kern$masklen)])
    }
  }
  return (k)  
}



selprojXselprojKernCompute <- function (kern1, kern2, x, x2) {
  funcName <- paste(kern1$comp[[1]]$type, "X", kern2$comp[[1]]$type, "KernCompute", sep="")
  func <- get(funcName, mode="function")

  mask <- .selprojKernMaskCombine(kern1$expmask, kern2$expmask)
  
  if ( nargs()>3 ) {
    d1 <- dim(as.array(x))[1]
    d2 <- dim(as.array(x2))[1]
    k <- array(0, dim=c(d1, d2))

    if (dim(as.array(x))[2] > 1 && dim(as.array(x2))[2] > 1 &&
        !any(is.na(mask))) {
      I1 <- .selprojDataMaskCompare(mask, x[,seq(kern1$masklen)])
      I2 <- .selprojDataMaskCompare(mask, x2[,seq(kern1$masklen)])

      if (!any(I1) || !any(I2))
        return (k)

      k[I1,I2] <- func(kern1$comp[[1]], kern2$comp[[1]], x[I1,-seq(kern1$masklen)], x2[I2,-seq(kern2$masklen)])
    }
  } else {
    d1 <- dim(as.array(x))[1]
    k <- array(0, dim=c(d1, d1))

    if (dim(as.array(x))[2] > 1 && !any(is.na(mask))) {
      I1 <- .selprojDataMaskCompare(mask, x[,seq(kern1$masklen)])

      if (!any(I1))
        return (k)

      k[I1,I1] <- func(kern1$comp[[1]], kern2$comp[[1]], x[I1,-seq(kern1$masklen)])
    }
  }
  return (k)  
}



selprojKernGradient <- function (kern, x, x2, covGrad) {
  funcName <- paste(kern$comp[[1]]$type, "KernGradient", sep="")
  func <- get(funcName, mode="function")

  if ( nargs()<4 ) {
    covGrad <- x2

    d1 <- dim(as.array(x))[1]

    stopifnot(dim(as.array(x))[2] > 1)

    I1 <- .selprojDataMaskCompare(kern$expmask, x[,seq(kern$masklen)])

    if (!any(I1))
      return (array(0, dim=kern$nParams))

    g <- func(kern$comp[[1]], x[I1,-seq(kern$masklen)], covGrad[I1,I1])
  } else {
    d1 <- dim(as.array(x))[1]
    d2 <- dim(as.array(x2))[1]

    stopifnot(dim(as.array(x))[2] > 1 && dim(as.array(x2))[2] > 1)
    I1 <- .selprojDataMaskCompare(kern$expmask, x[,seq(kern$masklen)])
    I2 <- .selprojDataMaskCompare(kern$expmask, x2[,seq(kern$masklen)])

    if (!any(I1) || !any(I2))
      return (array(0, dim=kern$nParams))

    g <- func(kern$comp[[1]], x[I1,-seq(kern$masklen)], x2[I2,-seq(kern$masklen)], covGrad[I1,I2])
  }

  return (g)
}


selprojXselprojKernGradient <- function (kern1, kern2, x, x2, covGrad) {
  funcName <- paste(kern1$comp[[1]]$type, "X", kern2$comp[[1]]$type, "KernGradient", sep="")
  func <- get(funcName, mode="function")
  mask <- .selprojKernMaskCombine(kern1$expmask, kern2$expmask)

  if ( nargs()<5 ) {
    covGrad <- x2

    d1 <- dim(as.array(x))[1]

    if (dim(as.array(x))[2] > 1 && !any(is.na(mask))) {
      I1 <- .selprojDataMaskCompare(mask, x[,seq(kern1$masklen)])

      if (!any(I1))
        return (list(g1=array(0, dim=kern1$nParams), g2=array(0, dim=kern2$nParams)))

      return (func(kern1$comp[[1]], kern2$comp[[1]], x[I1,-seq(kern1$masklen)], covGrad[I1,I1]))
    } else {
      return (list(g1=array(0, dim=kern1$nParams), g2=array(0, dim=kern2$nParams)))
    }
  } else {
    d1 <- dim(as.array(x))[1]
    d2 <- dim(as.array(x2))[1]

    if (dim(as.array(x))[2] > 1 && dim(as.array(x2))[2] > 1 &&
        !any(is.na(mask))) {
      I1 <- .selprojDataMaskCompare(mask, x[,seq(kern1$masklen)])
      I2 <- .selprojDataMaskCompare(mask, x2[,seq(kern1$masklen)])

      if (!any(I1) || !any(I2))
        return (list(g1=array(0, dim=kern1$nParams), g2=array(0, dim=kern2$nParams)))

      return (func(kern1$comp[[1]], kern2$comp[[1]], x[I1,-seq(kern1$masklen)], x2[I2,-seq(kern2$masklen)], covGrad[I1,I2]))
    } else {
      return (list(g1=array(0, dim=kern1$nParams), g2=array(0, dim=kern2$nParams)))
    }
  }
}


selprojKernDisplay <- function (kern, spaceNum=0) {
  spacing = matrix("", spaceNum+1)
  cat(spacing)
  cat("Selective projection kernel:\n")
  cat(spacing)
  cat("Experiment mask:", kern$expmask, "\n")
  for(i in seq(along=kern$comp)) 
    kernDisplay(kern$comp[[i]], spaceNum+2)
}
