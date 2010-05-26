modelGradientCheck <- function(model, ...) {

  params = modelExtractParam(model, only.values=FALSE)
  names = names(params)
  if (length(names) == 0) {
    for (i in 1:model$nParams)
      names[[i]] = paste('Param ', i, sep="")
  }

  if (length(names) != length(params))
    stop('Names array does not match length of params array')

  L = 0
  change = 1e-6
  origParams = params
  diff = matrix(0, 1, length(params))
  for (i in 1:length(params)) {
    params[i] = origParams[i] + change
    Lplus = modelObjective(params, model, ...)
    params[i] = origParams[i] - change
    Lminus = modelObjective(params, model, ...)
    diff[i] = (Lplus - Lminus)/(2*change)
    params[i] = origParams[i]
  }

  anal = modelGradient(origParams, model, ...)
  delta = anal-diff

  paramMaxDiff = max(abs(diff-anal))
  if (paramMaxDiff > 100*change) {
    l = 0
    for (i in 1:length(names)) {
      if (l < length(names[[i]]))
	l = length(names[[i]])
    }
    
    print('analytic   diffs     delta')
    for (i in 1:length(names)) {
      if(abs(delta[i]/max(c(abs(anal[i]), 1))) >= 1e-4) {
	spaceLen = l - length(names[[i]])
# 	space = char(repmat(32, 1, spaceLen));
	print(paste(names[[i]],':', round(c(anal[i], diff[i], diff[i]-anal[i]),6), sep=""))
      }
    }
  }
  print('Param max diff: ', round(paramMaxDiff,6))
}
