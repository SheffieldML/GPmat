\name{scoreList-class}
\Rdversion{1.1}
\docType{class}
\alias{scoreList-class}
\alias{show,scoreList-method}

\title{Class "scoreList"}
\description{  'scoreList' is an object which contain the genes, parameters,
  log-likelihoods and arguments of models. With the data in a scoreList
  item and the original data used for creating the models, the models
  can be reconstructed with the function 'generateModels'.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{scoreList(params, LLs, genes, modelArgs, knownTargets, TF)}.
}
\section{Slots}{
  \describe{
    \item{\code{params}:}{The parameters of the models.}
    \item{\code{LLs}:}{The log-likelihoods of the models.}
    \item{\code{genes}:}{The genes used in the models.}
    \item{\code{modelArgs}:}{A list of arguments used to generate the models.}
    \item{\code{knownTargets}:}{The list of known targets used in the ranking.}
    \item{\code{TF}:}{The TF used in the ranking.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "scoreList")}: Print a summary
      of the list}
  }
}
%\references{ ~put references to the literature/web site here ~ }
\author{Jonatan Ropponen, Antti Honkela}
%\note{ ~~further notes~~ }
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%
\seealso{
  \code{\link{GPrankTargets}, \link{GPrankTFs}, \link{generateModels}}.
}
%\seealso{
%	~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%	or \code{\linkS4class{CLASSNAME}} for links to other classes
%}
\examples{
showClass("scoreList")
}
\keyword{classes}
