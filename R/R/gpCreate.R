gpCreate <- function(q, d, X, y, options) {
  
  if (dim(X)[2]!=q)
    stop("Input matrix X does not have dimension ",q)
  if (dim(y)[2]!=d)
    stop("Input matrix y does not have dimension ",d)

  if (any(is.nan(y)) && !options$isMissingData)
    stop("NaN values in y, but no missing data declared.")
  if (options$isMissingData && options$isSpherical)
    stop("If there is missing data, spherical flag cannot be set.")

  y = as.matrix(y); X = as.matrix(X)

  model <- list(type="gp", y=y, X=X, approx=options$approx, beta = options$beta,
    learnScales=options$learnScales, scaleTransform=optimiDefaultConstraint("positive"),
    optimiseBeta=options$optimiseBeta, betaTransform=optimiDefaultConstraint("positive"),
    q=dim(X)[2], d=dim(y)[2], N=dim(y)[1])

#   ## Set up a mean function if one is given.
#   if (("meanFunction" %in% names(options)) && length(options$meanFunction)>0) {
#     if (is.list(options$meanFunction))
#       model$meanFunction = options$meanFunction
#     else
#       model$meanFunction = modelCreate(options$meanFunction, model$q, model$d,
# 			      options$meanFunctionOptions)
#   }

  model$optimiser = options$optimiser
  model$isSpherical = options$isSpherical
  model$isMissingData = options$isMissingData

  model$scale = matrix(1, 1, model$d)
  if (!model$isMissingData) {
    model$bias = colMeans(y)
  } else {
    for (i in 1:model$d) {
      model$indexPresent[[i]] = which(!is.nan(y[,i]))
      if (length(model$indexPresent[[i]])==0) {
	model$bias[i] = 0
# 	model$scale[i] = 1
      } else {
	model$bias[i] = mean(model$y[model$indexPresent[[i]], i])
# 	model$scale[i] = 1
      }
    }
  }
  
  if (("scale2var1" %in% names(options)) && (options$scale2var1)) {
    model$scale = sd(model$y)
    model$scale[which(model$scale == 0)] = 1
    if (model$learnScales)
      warning("Both learn scales and scale2var1 set for GP")
    if ("scaleVal" %in% names(options))
      warning("Both scale2var1 and scaleVal set for GP")
  }
  
  if("scaleVal" %in% names(options))
    model$scale = kronecker(matrix(1,1,model$d), options$scaleVal)
    # repmat(options$scaleVal, 1, model$d)

  model$m = gpComputeM(model)
  model$computeS = FALSE
  if (options$computeS) {
    model$computeS = TRUE
    model$S = model$m %*% t(model$m)
    if (model$approx != "ftc") # !strcmp(model.approx, 'ftc')
      stop("If compute S is set, approximation type must be 'ftc'.")
  }

  ## nParams is a sure way to tell if the kernel structure has been initialized
  if (is.list(options$kern) && ("nParams" %in% options$kern))
    model$kern = options$kern
  else
# browser()
    model$kern = kernCreate(model$X, options$kern)


#   if ("noise" %in% names(options)) {
#     if (is.list(options$noise)) #isstruct(options.noise)
#       model$noise = options$noise
#     else
#       model$noise = noiseCreate(options$noise, y)
#     
#     ## Set up noise model gradient storage.
#     model$nu = matrix(0, dim(y)[1], dim(y)[2])
#     model$g = matrix(0, dim(y)[1], dim(y)[2])
#     model$gamma = matrix(0, dim(y)[1], dim(y)[2])
#     
#     ## Initate noise model
#     model$noise = noiseCreate(noiseType, y) ## bug: noiseType has no value
#     
#     ## Set up storage for the expectations
#     model$expectations$f = model$y
#     model$expectations$ff = matrix(1, dim(model$y)[1], dim(model$y)[2])
#     model$expectations$fBar = matrix(1, dim(model$y)[1], dim(model$y)[2])
#     ## bug: numData has no value
#     model$expectations$fBarfBar = array(1,dim=(c(numData, numData, dim(model$y)[2])))
#   }


  if (options$approx == "ftc") {
    model$k = 0
    model$X_u = list()
    if (model$optimiseBeta && length(options$beta)==0) 
	stop("options.beta cannot be empty if it is being optimised.")
  } else if (options$approx %in% c("dtc", "dtcvar", "fitc", "pitc")) {
    ## Sub-sample inducing variables.
    model$k = options$numActive
    model$fixInducing = options$fixInducing
    if (options$fixInducing) {
      if (length(options$fixIndices) != options$numActive) {
	stop(paste("Length of indices for fixed inducing variables must ",
	  "match number of inducing variables"))
      }
      model$X_u = model$X[options$fixIndices, ]
      model$inducingIndices = options$fixIndices
    } else {
      ind = sample(1:model$N, size=model$N) #randperm(model$N)
      ind = ind[1:model$k]
      model$X_u = model$X[ind, ,drop=FALSE]
    }
  }

  if (model$k > model$N)
    stop("Number of active points cannot be greater than number of data.")

  if (model$approx == "pitc") { #strcmp(model.approx, 'pitc')
    numBlocks = ceiling(model$N/model$k)
    numPerBlock = ceiling(model$N/numBlocks)
    startVal = 1
    endVal = model$k
    model$blockEnd = matrix(0, 1, numBlocks)
    for (i in 1:numBlocks) {
      model$blockEnd[i] = endVal
      endVal = numPerBlock + endVal
      if (endVal > model$N)
	endVal = model$N
    }
  }

  initParams = gpExtractParam(model)

  ## This forces kernel computation.
# browser()
  model = gpExpandParam(model, initParams)

  return (model)
}
