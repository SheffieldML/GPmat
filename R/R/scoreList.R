setClass("scoreList", 
	representation(params = "list", LLs = "array", genes = "list", useGPsim = "array"))

	#prototype = list(params = list(), LLs = array(), genes = list(), useGPsim = array()))


scoreList <- function(params = list(), LLs = array(), genes = list(), useGPsim = array()) {

  new("scoreList", params = params, LLs = LLs, genes = genes, useGPsim = useGPsim)
}


is.scoreList <- function(object) {

  if (!is.na(charmatch(class(object), "scoreList"))) {
    return (TRUE)
  }  
  else {
    return (FALSE)
  }
}
