\name{scoreList-class}
\Rdversion{1.1}
\docType{class}
\alias{scoreList-class}
\alias{[,scoreList-method}
\alias{c,scoreList-method}
\alias{genes,scoreList-method}
\alias{genes<-,scoreList-method}
\alias{knownTargets,scoreList-method}
\alias{knownTargets<-,scoreList-method}
\alias{length,scoreList-method}
\alias{loglikelihoods,scoreList-method}
\alias{loglikelihoods<-,scoreList-method}
\alias{modelArgs,scoreList-method}
\alias{modelArgs<-,scoreList-method}
\alias{params,scoreList-method}
\alias{params<-,scoreList-method}
\alias{sharedModel,scoreList-method}
\alias{sharedModel<-,scoreList-method}
\alias{show,scoreList-method}
\alias{sort,scoreList-method}
\alias{TF,scoreList-method}
\alias{TF<-,scoreList-method}

\title{Class "scoreList"}
\description{  'scoreList' is an object which contain the genes, parameters,
  log-likelihoods and arguments of models. With the data in a scoreList
  item and the original data used for creating the models, the models
  can be reconstructed with the function 'generateModels'.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{scoreList(params, LLs,
  genes, modelArgs, knownTargets, TF, sharedModel)}.
}
\section{Slots}{
  \describe{
    \item{\code{params}:}{The parameters of the models.}
    \item{\code{LLs}:}{The log-likelihoods of the models.}
    \item{\code{genes}:}{The genes used in the models.}
    \item{\code{modelArgs}:}{A list of arguments used to generate the models.}
    \item{\code{knownTargets}:}{The list of known targets used in the ranking.}
    \item{\code{TF}:}{The TF used in the ranking.}
    \item{\code{sharedModel}:}{Shared model for known targets.}
  }
}
\section{Methods}{
  Class-specific methods:
  \describe{
    \item{\code{genes(object)}, \code{genes(object)<- value}}{Access and set \code{genes}}
    \item{\code{knownTargets(object)}, \code{knownTargets(object)<- value}}{Access and set \code{knownTargets}}
    \item{\code{loglikelihoods(object)}, \code{loglikelihoods(object)<- value}}{Access and set \code{loglikelihoods}}
    \item{\code{modelArgs(object)}, \code{modelArgs(object)<- value}}{Access and set \code{modelArgs}}
    \item{\code{params(object)}, \code{params(object)<- value}}{Access and set \code{params}}
    \item{\code{sharedModel(object)}, \code{sharedModel(object)<- value}}{Access and set \code{sharedModel}}
    \item{\code{TF(object)}, \code{TF(object)<- value}}{Access and set \code{TF}}
  }

  Standard generic methods:
  \describe{
    \item{\code{object[(index)}}{Conducts subsetting of the scoreList.}
    \item{\code{c(object, ...)}}{Concatenates scoreLists.}
    \item{\code{length(object)}}{Returns the length of the list.}
    \item{\code{show(object)}}{Informatively display object contents.}
    \item{\code{sort(object, decreasing=FALSE)}}{Sort the list according
  to log-likelihood}
  }
}
%\references{ ~put references to the literature/web site here ~ }
\author{Antti Honkela, Jonatan Ropponen}
%\note{ ~~further notes~~ }
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%
\seealso{
  \code{\link{GPRankTargets}}, \code{\link{GPRankTFs}},
  \code{\link{generateModels}}.
}
%\seealso{
%	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%	or \code{\linkS4class{CLASSNAME}} for links to other classes
%}
\examples{
showClass("scoreList")
}
\keyword{classes}
