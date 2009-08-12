printScoreListGenes <- function(scoreList) {

  for (i in 1:length(scoreList$data)) {
    genes <- scoreList$data[i][[1]]$genes
    LL <- scoreList$LLs[i]
    text <- c(genes, LL)
    cat(text)
    cat("\n")
  }
}

writeScoreListGenes <- function(scoreList, fileName) {

  text <- ""
  for (i in 1:length(scoreList$data)) {
    genes <- scoreList$data[i][[1]]$genes
    LL <- scoreList$LLs[i]
    newText <- c(genes, LL, "\n")
    text <- c(text, newText)
  }
  write(text, file = fileName)
}

writeGenes <- function(genes, fileName) {

  text <- ""
  for (i in 1:length(genes)) {
    text <- c(text, genes[i])
  }
  write(text, file = fileName)
}
