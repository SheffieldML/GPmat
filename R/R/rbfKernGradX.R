rbfKernGradX <- function (kern, x1, x2) {

  gX = array(0, c(dim(as.array(x2))[1], dim(as.array(x2))[2], dim(as.array(x1))[1]))
  for (i in 1:dim(x1)[1]) {
    gX[, , i] = rbfKernGradXpoint(kern, x1[i, ], x2)
  }

  return (gX)
}

rbfKernGradXpoint <- function (kern, x1, x2) {
  if (is.vector(x1))
    x1 = t(x1)

  gX = matrix(0, dim(as.array(x2))[1], dim(as.array(x2))[2])
  n2 = .dist2(x2, x1)
  wi2 = 0.5 * kern$inverseWidth
  rbfPart = kern$variance * exp(-n2 * wi2)
  for (i in 1:dim(x1)[2]) {
    gX[, i] = kern$inverseWidth * (x2[, i] - x1[i]) * rbfPart
  }

  if ('isNormalised' %in% names(kern) && (kern$isNormalised)) {
      gX = gX * sqrt(kern$inverseWidth/(2*pi))
  }

  return (gX)
}
