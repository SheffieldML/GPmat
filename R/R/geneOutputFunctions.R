printGenes <- function(scoreList) {

  for (i in 1:length(scoreList$data)) {
    genes <- scoreList$data[i][[1]]$genes
    LL <- scoreList$LLs[i]
    text <- c(genes, LL)
    cat(text)
    cat("\n")
  }
}

writeGenes <- function(scoreList, fileName) {

  text <- ""
  for (i in 1:length(scoreList$data)) {
    genes <- scoreList$data[i][[1]]$genes
    LL <- scoreList$LLs[i]
    newText <- c(genes, LL, "\n")
    text <- c(text, newText)
  }
  write(text, file = fileName)
}
