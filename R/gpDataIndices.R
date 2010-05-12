gpDataIndices <- function(model, dimNo, blockNo) {
  if (nargs() > 2) {
    if (model$approx != "pitc")
      stop("Block number only relevant for pitc approximation.")
    else {
      if (blockNo == 1)
	startVal = 1
      else
	startVal = model$blockEnd[blockNo-1]+1
      endVal = model$blockEnd[blockNo]
      if (model$isMissingData) {
	st = min(which(model$indexPresent[[dimNo]] >= startVal))
	fi = max(which(model$indexPresent[[dimNo]] <= endVal))
	ind = t(model$indexPresent[[dimNo]][st:fi])
      } else
	ind = startVal:endVal
    }
  } else {
    if (model$approx == "pitc")
      stop("Must give block number with PITC approximation.")
    else {
      if (model$isMissingData)
	ind = t(model$indexPresent[[dimNo]])
      else
	ind = 1:model$N
    }
  }

  return(ind)
}