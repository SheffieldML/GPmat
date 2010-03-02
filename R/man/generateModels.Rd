\name{generateModels}
\Rdversion{1.0}
\alias{generateModels}
\title{Generating models with the given data}
\description{
  'generateModels' recreates models based on the parameters
  stored in a scoreList.
}
\usage{
models <- generateModels(preprocData, scores)
}
\arguments{
  \item{preprocData}{The preprocessed data to be used.}
  \item{scores}{A scoreList object containing data of the models to be
    generated.}
}
%\details{
%  'generateModels' first creates a model with 'gpsimCreate' or 'gpdisimCreate'. Then 'gpsimExpandParam' or 'gpdisimExpandParam' is used to generate the original models without having to optimize them again.
%}
\value{
  'generateModels' returns a list of the generated models.
}
\author{
  Antti Honkela
}
\seealso{
  \code{\link{GPLearn}, \link{GPRankTargets}, \link{GPRankTFs}, \linkS4class{scoreList}}.
}
\examples{
## missing
}
\keyword{model}
