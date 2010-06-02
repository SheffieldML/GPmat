gpTest <- function(q=2, d=3, N=10, Nseq=4, k=5) {
  
  modelRet = list()
  kernType = list(type='cmpnd',comp=list('rbf', 'lin', 'rbfard', 'mlp', 'mlpard', 'white'))
  kernType = list(type='cmpnd',comp=list('rbf', 'white'))
  meanFunctionType = 'mlp'
  learnScales = TRUE ## test learning of output scales.
#   learnScales = FALSE
  X = matrix(rnorm(N*q), N, q)
  Yorig = matrix(rnorm(N*d), N, d)
  indMissing = which(matrix(rnorm(N*q), N, q) > 0.7)
  approxType = list("ftc", "dtc", "dtcvar", "fitc", "pitc")
#   approxType = list("pitc")
  counter = 0

  for (optimiseBeta in FALSE:TRUE) {
#   for (optimiseBeta in FALSE) {

#   for (meanFunction in FALSE:TRUE) {
    for (meanFunction in FALSE) {

      for (missing in FALSE:TRUE) {
#       for (missing in FALSE) {

	for (fixInducing in FALSE:TRUE) {
# 	for (fixInducing in FALSE) {

	  Y = Yorig

	  if (missing)
	    Y[indMissing] = NaN

	  if (meanFunction && missing)
	    next

	  for (a in 1:length(approxType)) {
	    options = gpOptions(approxType[[a]])
	    options$learnScales = learnScales
	    options$kern = kernType
	    options$numActive = k
	    options$isSpherical = !missing
	    options$isMissingData = missing
	    options$fixInducing = fixInducing
	    options$optimiseBeta = optimiseBeta
	    if (optimiseBeta && approxType[[a]]=="ftc")
	      options$beta = 1000
	    else if (missing && approxType[[a]]=="dtcvar")
	      next

	    if (meanFunction)
	      print(paste("Mean Function installed, with ", approxType[[a]],
		" approximation.",sep=""))
	    else
	      print(paste(approxType[[a]], " approximation.",sep=""))

	    if (missing) {
	      print("Missing data used.")
	    }

	    if (fixInducing) {
	      print("Inducing variables fixed.")
	      options$fixIndices = round(seq(1, dim(Y)[1], len=k))
	    }
	    if (!optimiseBeta)
	      print("Beta not optimised.")
	  
	    if (meanFunction) {
	      options$meanFunction = meanFunctionType
	      options$meanFunctionOptions =
		get(paste(meanFunctionType, "Options",sep=""), mode="function")()
	    }

	    model = gpCreate(q, d, X, Y, options)

	    initParams = gpExtractParam(model)
	    ## this creates some nasty parameters.
	    initParams=matrix(rnorm(prod(dim(as.matrix(initParams)))),
			dim(as.matrix(initParams))[1], dim(as.matrix(initParams))[2])
	    # / (100*matrix(rnorm(prod(dim(as.matrix(initParams)))),
	    # 		dim(as.matrix(initParams))[1], dim(as.matrix(initParams))[2])

	    ## This forces kernel computation.
	    model = gpExpandParam(model, initParams)
	    gpCovGradsTest(model)
	    modelGradientCheck(model)
	    counter = counter + 1
	    modelRet[[counter]] = model
	  }
	}
      }
    }
  }
  return (modelRet)
}