gpPosteriorMeanVar <- function(model, X, varsigma.return=FALSE) {

  if (!"alpha" %in% names(model))
    model = gpComputeAlpha(model)

  maxMemory = 1000000

  if (model$approx == "ftc")
    chunkSize = ceiling(maxMemory/model$N)
  else if (model$approx %in% c("dtc", "dtcvar", "fitc", "pitc"))
    chunkSize = ceiling(maxMemory/model$k)

  mu = matrix(0, dim(X)[1], model$d)
#   if (varsigma.return) #if (nargout > 1)
#     varsigma = matrix(0, dim(X)[1], model$d)

  startVal = 1
  endVal = chunkSize
  if (endVal > dim(X)[1])
    endVal = dim(X)[1]

  while (startVal <= dim(X)[1]) {
    indices = startVal:endVal

    ## Compute kernel for new point.
    if (model$approx == "ftc")
      KX_star = kernCompute(model$kern, model$X, X[indices,]) #in kernFunctions.R
    else if (model$approx %in% c("dtc", "dtcvar", "fitc", "pitc"))
      KX_star = kernCompute(model$kern, model$X_u, X[indices,]) #in kernFunctions.R
    
    ## Compute mean, using precomputed alpha vector.
    if ((!"isMissingData" %in% names(model)) || !model$isMissingData || model$approx != "ftc")
      mu[indices,] = KX_star %*% model$alpha
    else {
      for (i in 1:model$d)
	mu[indices, i] = KX_star[model$indexPresent[[i]],]
	    %*% model$alpha[model$indexPresent[[i]], i]
    }
    
    ## Compute variances if requried.
    if (varsigma.return) { #if (nargout > 1)
      varsigma = matrix(0, dim(X)[1], model$d)
      if (!("isSpherical" %in% names(model)) || model$isSpherical) {
	## Compute diagonal of kernel for new point.
	diagK = kernDiagCompute(model$kern, X[indices,])
	if (model$approx == "ftc")
	  Kinvk = model$invK_uu %*% KX_star
	else if (model$approx %in% c("dtc", "dtcvar", "fitc", "pitc")) {
	  Kinvk = (model$invK_uu - (1/model$beta)%*%model$Ainv) %*% KX_star
	}
	varsig = diagK - t(colSums(KX_star * Kinvk, 1))
	if ("beta"  %in% names(model)) {
	  varsig = varsig + (1/model$beta)
	}
	varsigma[indices,] = kronecker(matrix(1,1,model$d), varsig)
      } else {
	diagK = kernDiagCompute(model$kern, X[indices,])
	for (i in 1:model$d) {
	  ind = model$indexPresent[[i]]
	  if (model$approx == "ftc")
	    Kinvk = model$invK_uu[[i]] %*% KX_star[ind,]
	  else {
	    stop(c("Non-spherical not yet implemented for any approximation",
		  "other than 'ftc'."))
	  }
	  varsigma[indices, i] = diagK - t(colSums(KX_star[ind,] * Kinvk, 1))
	}
      }
    }
    
    
    ## Rescale the mean
    mu[indices,] = mu[indices,] * kronecker(matrix(1,length(indices),1), model$scale)
    ## Add the bias back in.
    mu[indices,] = mu[indices,] + kronecker(matrix(1,length(indices),1), model$bias)
    ## If the mean function is present, add it it.
    if (("meanFunction" %in% names(model)) && length(model$meanFunction)>0) {
      mu[indices,] = mu[indices,] + modelOut(model$meanFunction, X[indices,])
    }
    ## rescale the variances
    if (varsigma.return) { #if (nargout > 1)
      varsigma[indices,] = varsigma[indices,] *
	kronecker(matrix(1,length(indices),1), model$scale*model$scale)
    }

    ## Prepare for the next chunk.
    startVal = endVal + 1
    endVal = endVal + chunkSize
    if (endVal > dim(X)[1])
      endVal = dim(X)[1]
  } #while

  if (varsigma.return)
    return (list(mu=mu, varsigma=varsigma))
  else
    return (mu)
}