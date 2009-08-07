load("preprocData.RData")
regulators <- readRegulators("possibleRegulatorsOfGata3.txt")
scoreList <- GPscoreListFixedTargets(preprocData, targets = "100924_at", searchedGenes = c("100924_at", regulators), search = TRUE)
save(scoreList, file = "~/R/GPanalysis/GPsimAnalysis1.1/Gata3TFsList.RData")
