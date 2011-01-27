\name{kernCreate}
\Rdversion{1.0}
\alias{kernCreate}
\alias{kernParamInit}
\alias{cmpndKernParamInit}
\alias{disimKernParamInit}
\alias{mlpKernParamInit}
\alias{multiKernParamInit}
\alias{rbfKernParamInit}
\alias{simKernParamInit}
\alias{translateKernParamInit}
\alias{whiteKernParamInit}
\alias{priorCreate}
\alias{gammaPriorParamInit}
\alias{invgammaPriorParamInit}
\alias{priorParamInit}

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

  The *KernParamInit functions perform initialisation specific to
  different types of kernels. They should not be called directly.
}
\value{
  \item{kern}{The kernel structure.}
}
\seealso{
\code{\link{kernDisplay}, \link{modelTieParam}}.
}
\examples{
# Create a multi kernel with two rbf blocks with bounded inverse widths
invWidthBounds <- c(0.5, 2)
kernType <- list(type="multi", comp=list())
for (i in 1:2)
  kernType$comp[[i]] <- list(type="parametric", realType="rbf",
                              options=list(isNormalised=TRUE,
                                inverseWidthBounds=invWidthBounds))
kern <- kernCreate(1, kernType)

# Tie the inverse with parameters of the component RBF kernels
kern <- modelTieParam(kern, list(tieWidth="inverseWidth"))
kernDisplay(kern)
}
\keyword{model}
