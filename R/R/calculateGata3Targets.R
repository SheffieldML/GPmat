library("gpsim", lib.loc = "~/R/Rlibs")
load("~/R/GPanalysis/GPsimAnalysis1.2/GPsimAnalysis1.2Package/R/preprocData.RData")
scoreList <- GPscoreListFixedTF(preprocData, TF = "100924_at")
save(scoreList, file = "~/R/GPanalysis/GPsimAnalysis1.2/GPsimAnalysis1.2Package/R/Gata3Targets.RData")
