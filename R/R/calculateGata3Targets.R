load("preprocData.RData")
scoreList <- GPscoreListFixedTF(preprocData, TF = "100924_at")
save(scoreList, file = "~/R/GPanalysis/GPsimAnalysis1.0/Gata3Targets.RData")
