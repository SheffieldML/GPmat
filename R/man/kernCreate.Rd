\name{kernCreate}
\Rdversion{1.0}
\alias{kernCreate}
\title{Initialise a kernel structure.}
\description{
  Initialise a kernel structure.
}
\usage{
kern <- kernCreate(X, type)
kern <- kernCreate(dim, type)
}
\arguments{
  \item{X}{Input data values (from which kernel will later be computed).}
  
  \item{type}{Type of kernel to be created, some standard types are
    'rbf', 'white', 'sim' and 'disim'. If a list of the form
    \code{list(type='cmpnd', comp=c('rbf', 'rbf', 'white'))} is used a
    compound kernel based on
    the sum of the individual kernels will be created.
    Parameters can be passed to kernels using type
    \code{list(type='parametric', options=list(opt=val), realType=...)},
    where \code{realType} is the type that would be used otherwise.
    }
  
  \item{dim}{input dimension of the design matrix (i.e. number of
    features in the design matrix).}
}
\details{
  \code{kern <- kernCreate(X, type)}
  input points and a kernel type.

  \code{kern <- kernCreate(dim, type)}
  creates a kernel matrix structure given the dimensions of the design matrix and the kernel type.

}
\value{
  \item{kern}{The kernel structure.}
}
\seealso{
\code{\link{kernParamInit}}.
}
\examples{
## missing
}
\keyword{model}
