\name{GPModel-class}
\Rdversion{1.1}
\docType{class}
\alias{GPModel-class}
\alias{is.GPModel}
\alias{modelStruct}
\alias{modelStruct<-}
\alias{modelType}
\alias{initialize,GPModel-method}
\alias{is.GPModel,GPModel-method}
\alias{modelStruct,GPModel-method}
\alias{modelStruct<-,GPModel,list-method}
\alias{modelType,GPModel-method}
\alias{show,GPModel-method}

\title{A container for gpsim models}
\description{The class is a container for the internal representation
  of models used by the \code{gpsim} package.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GPModel", model)}.
}
\section{Slots}{
  \describe{
    \item{\code{model}:}{A model object used internally by the code of
      the \code{gpsim} package}
    \item{\code{type}:}{Type of the model object}
  }
}
\section{Methods}{
  \describe{
    \item{\code{modelStruct(object)}, \code{modelStruct(object)<-
	value}}{Access and set the internal \code{model} structure}
    \item{\code{modelType(object)}}{Access the internal \code{type} values}
    \item{\code{show(object)}}{Informatively display object contents.}
    \item{\code{is.GPModel(object)}}{Check if object is a GPModel.}
    \item{\code{initialize("GPModel")}}{Object instantiation, used by new; not to be called directly by the user.}
  }
}
%\references{ ~put references to the literature/web site here ~ }
\author{Antti Honkela, Jonatan Ropponen}
%\note{ ~~further notes~~ }
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%
\seealso{
  \code{\link{GPLearn}}, \code{\link{GPRankTargets}}, \code{\link{GPRankTFs}},
  \code{\link{generateModels}}, \code{\link{modelExtractParam}},
  \code{\link{modelLogLikelihood}}.
}
\examples{
showClass("GPModel")
}
\keyword{classes}
