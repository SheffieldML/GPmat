multiKernFixBlocks <- function(kern, blocks = 1) {

  source("multiKernCacheBlock.R")

  length <- 30
  op <- options(digits.secs=15)
  time <- Sys.time()
  uuid <- ""
  for (i in 1:length) {
    char <- substr(time, i, i)
    if (is.na(charmatch(char, ".")) && is.na(charmatch(char, ":")) && is.na(charmatch(char, "-")) && is.na(charmatch(char, " "))) {
      uuid <- paste(uuid, char, sep="")
    }
  }
  uuid <- paste('a', uuid, sep="")

  kern$fixedBlocks <- array(0, dim = c(1, kern$numBlocks))
  kern$fixedBlocks[blocks] <- 1
  kern$uuid <- uuid
  multiKernCacheBlock(kern)

  return(kern)
}
