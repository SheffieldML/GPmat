\name{GPplot}
\Rdversion{1.0}
\alias{GPplot}
\title{Drawing plots of expression levels}
\description{
  Draws plots of the expression levels of specified genes.
}
\usage{
GPplot(data, savepath = '', doParams = FALSE, selectGenes = data$genes, nameMapping = NULL)
}
\arguments{
  \item{data}{The data used for drawing the figures. The data is received as output from the function GPrank. The data is assumed to have one trancription factor and one target in each model.}
  \item{savepath}{The location in the file system where the images are saved.}
  \item{doParams}{}
  \item{selectGenes}{The genes for which the figures are drawn. By default, all genes of the data are selected.}
  \item{nameMapping}{The function used for mapping the names of the genes for the figures.}
}
\details{
  The function draws 3 plots in each figure. One of the plots represents the expression level of the transcription factor, one the expression level of the target and one the protein concentration of the inferred transcription factor.
}
\author{
  Jonatan Ropponen
}
\seealso{
  \code{\link{GPrank}}.
}
\examples{
## missing
}
\keyword{model}
