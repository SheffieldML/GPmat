multiKernCacheBlock <- function(kern, fhandle = NULL, arg = NULL, i = NULL, j = NULL, X = NULL, X2 = NULL) {

#  oneArgument <- FALSE
#
#  if(is.null(fhandle) && is.null(arg) && is.null(i) && is.null(j) && is.null(X) && is.null(X2)) {
#    oneArgument <- TRUE
#  }

  if(!exists("cacheEnvir") || !exists("cache")) {
    cacheEnvir <- new.env("cacheEnvir")
    cache <- list()
    cache$kern$uuid <- array(dim = kern$numBlocks)
    assign("cache", cache, pos = cacheEnvir)
    K <- array()
    return (K)
  }

  cache <- get("cache", envir = as.environment(cacheEnvir)) 
  myCache <- cache$kern$uuid[i,j]
  key <- rbind(X, X2)

  for (k in 1:length(myCache)) {
    if (dim(key)[2] == length(myCache[k][1]) && all(key == myCache[k][1])) {
      K = myCache[k][2]
      return (K)
    }
  }

  # No match if we get all the way here
  if (is.null(X2)) {
    K = fhandle(arg, X)
  }
  else {
    K = fhandle(arg, X, X2)
  }
  cache$kern$uuid[i, j][end+1] <- c(key, K)
}
