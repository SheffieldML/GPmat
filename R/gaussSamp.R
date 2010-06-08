gaussSamp <- function(Sigma, numSamps) {

  eigVecs = eigen(Sigma)
  U = eigVecs$vectors; V = eigVecs$values
  dims = dim(Sigma)[1]
  y = matrix(rnorm(numSamps*dims), numSamps, dims) #randn(numSamps, dims)
  V[V<0] = as.complex(V[V<0])
  y = y %*% diag(sqrt(V)) %*% t(U)

  return (y)
}