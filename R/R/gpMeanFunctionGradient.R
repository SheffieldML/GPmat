gpMeanFunctionGradient <- function(model) {

  if ("isSpherical" %in% names(model) && !model$isSpherical)
    stop("Currently only implemented for spherical")
  else
    if (model$isMissingData)
      stop("Currently not implemented for missing data.")
    
  if ("meanFunction" %in% names(model) && length(model$meanFunction)>0) {
    g = matrix(0, 1, model$meanFunction$numParams)
    ## compute gradients here.
    if (model$approx == "ftc")
      gmu = model$invK_uu%*%model$m
    else if (model$approx %in% c("dtc", "dtcvar"))
      gmu = (model$m - t(model$K_uf)%*%model$Ainv%*%(model$K_uf%*%model$m))*model$beta
    else if (model$approx == "fitc") {
      Dinvm = model$Dinv%*%model$m
      gmu = (Dinvm-(model$Dinv %*% t(model$K_uf))
	      %*%(model$Ainv%*%model$K_uf)%*%Dinvm)*model$beta
    } else if (model$approx == "pitc") {
      ## Loop through the blocks computing each part to be added.
      gmu = matrix(0, model$N, model$d)
      K_ufDinvm = matrix(0,model$k, model$d)
      K_ufDinv = matrix(0,model$k, model$N)
      for (i in 1:length(model$blockEnd)) {
	ind = gpBlockIndices(model, i)
	Dinvm[[i]] = model$Dinv[[i]]%*%model$m[ind, ,drop=F]
	K_ufDinvm = K_ufDinvm + model$K_uf[, ind,drop=F]%*%Dinvm[[i]]
      }
      for (i in 1:length(model$blockEnd)) {
	ind = gpBlockIndices(model, i)
	gmu[ind, ] = (Dinvm[[i]] - model$Dinv[[i]]
		      %*%t(model$K_uf[, ind,drop=F])%*%(model$Ainv%*%K_ufDinvm))*model$beta
      }
    }
    
    gmu = gmu/kronecker(matrix(1,model$N, 1),model$scale)
    goutputDparam = modelOutputGrad(model$meanFunction, model$X)
    for (i in 1:model$meanFunction$numParams)
      g[1, i] = sum(gmu * drop(goutputDparam[, i, ])) # drop=squeeze
  } else
    g = list()

  return (g)
}