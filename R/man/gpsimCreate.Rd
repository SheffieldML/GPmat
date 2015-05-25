\name{gpsimCreate}
\Rdversion{1.0}
\alias{gpsimCreate}
\alias{gpdisimCreate}
\title{Create a GPSIM/GPDISIM model.}
\description{
  creates a model for single input motifs with Gaussian processes.
}
\usage{
  gpsimCreate(Ngenes, Ntf, times, y,
              yvar, options, genes=NULL, annotation=NULL)
  gpdisimCreate(Ngenes, Ntf, times, y,
                yvar, options, genes=NULL, annotation=NULL)
}
\arguments{
  \item{Ngenes}{number of genes to be modelled in the system.}
  \item{Ntf}{number of proteins to be modelled in the system.}
  \item{times}{the time points where the data is to be modelled.}
  \item{y}{the values of each gene at the different time points.}
  \item{yvar}{the variances of each gene at the different time points.}
  \item{options}{options structure (optional).}
  \item{genes}{names of the probes the model is for}
  \item{annotation}{(optional) annotation for the probe names}
}
\value{
  \item{model}{model structure containing default parameterisation.}
}
\details{
  These functions are meant to be used through \code{\link{GPLearn}}.
}
\seealso{
\code{\link{modelExtractParam}, \link{modelOptimise}, \link{GPLearn}}.
}
\examples{
## missing, see GPLearn
}
\keyword{model}
