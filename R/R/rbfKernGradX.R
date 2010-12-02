rbfKernGradX <- function (kern, X, X2) {

  gX = array(0, c(dim(as.array(X2))[1], dim(as.array(X2))[2], dim(as.array(X))[1]))
  for (i in 1:dim(X)[1]) {
    gX[, , i] = rbfKernGradXpoint(kern, X[i, ], X2)
  }

  return (gX)
}

rbfKernGradXpoint <- function (kern, x, X2) {
  if (is.vector(x))
    x = t(x)

  gX = matrix(0, dim(as.array(X2))[1], dim(as.array(X2))[2])
  n2 = .dist2(X2, x)
  wi2 = 0.5 * kern$inverseWidth
  rbfPart = kern$variance * exp(-n2 * wi2)
  for (i in 1:dim(x)[2]) {
    gX[, i] = kern$inverseWidth * (X2[, i] - x[i]) * rbfPart
  }

  if ('isNormalised' %in% names(kern) && (kern$isNormalised)) {
      gX = gX * sqrt(kern$inverseWidth/(2*pi))
  }

  return (gX)
}
