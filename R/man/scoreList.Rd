\name{scoreList}
\Rdversion{1.0}
\alias{scoreList}
\title{Containing data for reconstructing models}
\description{
  'scoreList' is a S4 class the items of which contain the genes, parameters, log-likelihoods and useGPsim values of models. With the data in a scoreList item and the original data used for creating the models, the models can be reconstructed with the function 'generateModels'.

  'is.scoreList' is a function determining whether the given object is a scoreList.
}
\usage{
list <- scoreList(params = list(), LLs = array(), genes = list(), useGPsim = array())
value <- is.scoreList(object)
}
\arguments{
  \item{params}{The parameters of the models.}
  \item{LLs}{The log-likelihoods of the models.}
  \item{genes}{The genes used in the models.}
  \item{useGPsim}{A logical value determining whether GPsim was used in forming the model; if false, GPdisim was used.}
  \item{object}{The object given to is.scoreList.}
}
\details{
  ScoreList items can be produced as output of the functions GPtankTargets and GPrankTFs.
}
\value{
  'is.scoreList' returns a logical value determining whether the given object is a scoreList.
}
\author{
  Jonatan Ropponen
}
\seealso{
  \code{\link{GPrank}, \link{formModel}, \link{GPrankTFs}, \link{GPrankTargets}, \link{generateModels}}.
}
\examples{
## missing
}
\keyword{model}
