\name{scoreList-class}
\Rdversion{1.1}
\docType{class}
\alias{scoreList-class}
\alias{write.scores}
\alias{genes}
\alias{genes<-}
\alias{knownTargets}
\alias{knownTargets<-}
\alias{loglikelihoods}
\alias{loglikelihoods<-}
\alias{baseloglikelihoods}
\alias{baseloglikelihoods<-}
\alias{modelArgs}
\alias{modelArgs<-}
\alias{params}
\alias{params<-}
\alias{sharedModel}
\alias{sharedModel<-}
\alias{TF}
\alias{TF<-}
\alias{[,scoreList-method}
\alias{c,scoreList-method}
\alias{write.scores,scoreList-method}
\alias{genes,scoreList-method}
\alias{genes<-,scoreList,list-method}
\alias{knownTargets,scoreList-method}
\alias{knownTargets<-,scoreList,character-method}
\alias{length,scoreList-method}
\alias{loglikelihoods,scoreList-method}
\alias{loglikelihoods<-,scoreList,numeric-method}
\alias{baseloglikelihoods,scoreList-method}
\alias{baseloglikelihoods<-,scoreList,numeric-method}
\alias{modelArgs,scoreList-method}
\alias{modelArgs<-,scoreList,list-method}
\alias{params,scoreList-method}
\alias{params<-,scoreList,list-method}
\alias{sharedModel,scoreList-method}
\alias{sharedModel<-,scoreList,list-method}
\alias{show,scoreList-method}
\alias{sort,scoreList-method}
\alias{TF,scoreList-method}
\alias{TF<-,scoreList,character-method}

\title{Class "scoreList"}
\description{  'scoreList' is an object which contain the genes, parameters,
  log-likelihoods and arguments of models. With the data in a scoreList
  item and the original data used for creating the models, the models
  can be reconstructed with the function 'generateModels'.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{scoreList(params, loglikelihoods,
  genes, modelArgs, knownTargets, TF, sharedModel)}.
}
\section{Slots}{
  \describe{
    \item{\code{params}:}{The parameters of the models.}
    \item{\code{loglikelihoods}:}{The log-likelihoods of the models.}
    \item{\code{baseloglikelihoods}:}{The log-likelihoods of
      corresponding null models.}
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
    \item{\code{write.scores(object, ...)}}{Writes the log-likelihoods
      and null log-likelihoods. Accepts any options \code{write.table} does.}
    \item{\code{genes(object)}, \code{genes(object)<- value}}{Access and set \code{genes}}
    \item{\code{knownTargets(object)}, \code{knownTargets(object)<- value}}{Access and set \code{knownTargets}}
    \item{\code{loglikelihoods(object)}, \code{loglikelihoods(object)<- value}}{Access and set \code{loglikelihoods}}
    \item{\code{baseloglikelihoods(object)}, \code{baseloglikelihoods(object)<- value}}{Access and set \code{baseloglikelihoods}}
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
  \code{\link{generateModels}}, \code{\link{write.table}}.
}
%\seealso{
%	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%	or \code{\linkS4class{CLASSNAME}} for links to other classes
%}
\examples{
showClass("scoreList")
}
\keyword{classes}
