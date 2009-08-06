calculateGata3LLs <- function() {
  preprocData <- processData(mmgmos_exprs_refseq)
  scoreList <- GPscoreListFixedTF(preprocData, TF = "100924_at")
}
