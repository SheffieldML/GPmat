gpBlockIndices <- function(model, blockNo) {
  if (model$approx != "pitc")
    stop("Block number only relevant for pitc approximation.")
  else {
    if (blockNo == 1)
      startVal = 1
    else
      startVal = model$blockEnd[blockNo-1]+1
    
    endVal = model$blockEnd[blockNo]
    ind = startVal:endVal
  }

  return(ind)
}