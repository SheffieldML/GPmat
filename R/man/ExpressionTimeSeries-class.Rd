\name{ExpressionTimeSeries-class}
\Rdversion{1.1}
\docType{class}
\alias{ExpressionTimeSeries-class}
\alias{var.exprs}
\alias{var.exprs<-}
\alias{initialize,ExpressionTimeSeries-method}
\alias{var.exprs,ExpressionTimeSeries-method}
\alias{var.exprs<-,ExpressionTimeSeries-method}

\title{Class to contain time series expression assays}
\description{Container for time series expression assays and
  experimental metadata. \code{ExpressionTimeSeries} class is derived
  from \code{\linkS4class{ExpressionSet}}, and requires fields
  \code{experiments} and \code{modeltime} in \code{phenoData}.}
\section{Extends}{
Directly extends class \code{\linkS4class{ExpressionSet}}.
}
\section{Objects from the Class}{
  \code{new("ExpressionTimeSeries")}

  \code{new("ExpressionTimeSeries",
    phenoData = new("AnnotatedDataFrame"),
    featureData = new("AnnotatedDataFrame"),
    experimentData = new("MIAME"),
    annotation = character(0),
    protocolData = phenoData[,integer(0)],
    exprs = new("matrix"),
    var.exprs = new("matrix"))
  }

  This creates an \code{ExpressionTimeSeries} with \code{assayData} implicitly
  created to contain \code{exprs} and \code{var.exprs}.

  \code{new("ExpressionTimeSeries",
    assayData = assayDataNew(exprs=new("matrix")),
    phenoData = new("AnnotatedDataFrame"),
    featureData = new("AnnotatedDataFrame"),
    experimentData = new("MIAME"),
    annotation = character(0),
    protocolData = phenoData[,integer(0)])
  }

  This creates an \code{ExpressionTimeSeries} with \code{assayData} provided
  explicitly. In this form, the only required named argument is
  \code{assayData}.

  \code{ExpressionTimeSeries} instances are usually created through
  \code{new("ExpressionTimeSeries", ...)}. Usually the arguments to \code{new}
  include \code{exprs} (a matrix of expression data, with features
  corresponding to rows and samples to columns), \code{var.exprs},
  \code{phenoData},
  \code{featureData}, \code{experimentData}, \code{annotation}, and
  \code{protocolData}. 
  \code{phenoData}, \code{featureData}, \code{experimentData},
  \code{annotation}, and \code{protocolData} can be missing, in which case
  they are assigned default values.
}
\section{Slots}{
  \describe{
    \item{\code{assayData}:}{Inherited from
      \code{\linkS4class{ExpressionSet}}.
      The models in \code{gpsim} package assume that \code{exprs}
      contains absolute (i.e. non-logarithmic) expression values.
      The member \code{var.exprs} may contain variances of the values.}
    \item{\code{phenoData}:}{Inherited from
      \code{\linkS4class{ExpressionSet}}.
      The following fields are required: \code{experiments} which
      contains integers from 1 to N with measurements from the same
      biological assay having the same number; \code{modeltime}
      which contains observation times in model units.
    }
    \item{\code{featureData}:}{Inherited from \code{\linkS4class{ExpressionSet}}.}
    \item{\code{experimentData}:}{Inherited from \code{\linkS4class{ExpressionSet}}.}
    \item{\code{annotation}:}{Inherited from \code{\linkS4class{ExpressionSet}}.}
    \item{\code{protocolData}:}{Inherited from \code{\linkS4class{ExpressionSet}}.}
    \item{\code{.__classVersion__}:}{Inherited from \code{\linkS4class{ExpressionSet}}.}
  }
}
\section{Methods}{
  See also methods for \code{\linkS4class{ExpressionSet}}.
  \describe{
    \item{\code{var.exprs(object)}, \code{var.exprs(object)<- value}}{Access and set \code{var.exprs}}
    \item{\code{initialize("ExpressionTimeSeries")}}{Object instantiation, used by new; not to be called directly by the user.}
  }
}
\author{Antti Honkela, Jonatan Ropponen}
\seealso{
  \code{\link{processData}}, \code{\link{processRawData}}.
}
\examples{
showClass("ExpressionTimeSeries")
}
\keyword{classes}
