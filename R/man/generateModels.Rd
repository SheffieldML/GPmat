\name{generateModels}
\Rdversion{1.0}
\alias{generateModels}
\alias{searchExpressionData}
\title{Generating models with the given data}
\description{
  'generateModels' creates models based on the given parameters without having to optimize the models again. It's possible to use a scoreList iteam and preprocessed data for this or alternatively give the data individually.

  'searchExpressionData' searches for the data of specified genes from the expression and variance data.
}
\usage{
listOfModels <- generateModels(preprocData = NULL, y = NULL, yvar = NULL, allGenes = NULL, times = NULL, scoreList = NULL, params = NULL, genes = NULL, useGPsim = NULL)
newData <- searchExpressionData(y, yvar, allGenes, times, searchedGenes)
}
\arguments{
  \item{preprocData}{The preprocessed data to be used.}
  \item{y}{The expression data if not given in preprocData.}
  \item{yvar}{The variance data if not given in preprocData..}
  \item{allGenes}{An array of all genes in the expression data if not given in preprocData.}
  \item{times}{The time points used if not given in preprocData.}
  \item{scoreList}{A scoreList item containing data of the models to be generated.}
  \item{params}{The parameters of the models if not given in scoreList.}
  \item{genes}{The genes used in the models if not given in scoreList.}
  \item{useGPsim}{An array of the useGPsim values of the models. Each value is a logical value determining whether GPsim was used; if false, GPdisim was used.}
  \item{searchedGenes}{The genes the data of which is searched.}
}
\details{
  'generateModels' first creates a model with 'gpsimCreate' or 'gpdisimCreate'. Then 'gpsimExpandParam' or 'gpdisimExpandParam' is used to generate the original models without having to optimize them again.
}
\value{
  'generateModels' returns a list of the generated models.

  'searchExpressionData' returns the data of the specified genes.
}
\author{
  Jonatan Ropponen
}
\seealso{
  \code{\link{GPrank}, \link{formModel}, \link{GPrankTFs}, \link{GPrankTargets}, \link{scoreList}}.
}
\examples{
## missing
}
\keyword{model}
