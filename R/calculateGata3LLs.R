calculateGata3LLs <- function() {

  source("preprocess_affy.R")
  source("processData.R")
  source("GPscoreListFixedTF.R")

  preprocData <- processData(mmgmos_exprs_refseq)
  scoreList <- GPscoreListFixedTF(preprocData, TF = "100924_at")
}
