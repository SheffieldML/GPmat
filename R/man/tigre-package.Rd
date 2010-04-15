\name{tigre-package}
\alias{tigre-package}
\alias{tigre}
\docType{package}
\title{
  tigre - Transcription factor Inference through Gaussian process
  Reconstruction of Expression
}
\description{
  This package implements the method of Gao et al. (2008) and Honkela et
  al. (2010) for Gaussian process modelling single input motif
  regulatory systems with time-series expression data. The method can be
  used to rank potential targets of transcription factors based on such
  data.
}
\details{
  \tabular{ll}{
    Package: \tab tigre\cr
    Type: \tab Package\cr
    Version: \tab 0.9.4\cr
    Date: \tab 2010-04-12\cr
    License: \tab A-GPL Version 3\cr
  }
  For details of using the package please refer to the Vignette.
}
\author{
  Antti Honkela, Pei Gao, Jonatan Ropponen, Magnus Rattray, Neil D. Lawrence

  Maintainer: Antti Honkela <antti.honkela@tkk.fi>
}
\references{
  P.~Gao, A.~Honkela, M.~Rattray, and N.~D.~Lawrence.
  Gaussian process modelling of latent chemical species:
  applications to inferring transcription factor activities.
  \emph{Bioinformatics} 24(16):i70--i75, 2008.

  A.~Honkela, C.~Girardot, E.~H. Gustafson, Y.-H. Liu,
  E.~E.~M. Furlong, N.~D. Lawrence, and M.~Rattray.
  Model-based method for transcription factor target
  identification with limited data.
  \emph{Proceedings of the National Academy of Sciences of the USA}
  (2010). In press.
}
\keyword{ package }
\seealso{
  \code{\link[puma:puma-package]{puma}}
}
\examples{\dontrun{
  # Load a mmgmos preprocessed fragment of the Drosophila developmental
  # time series
  data(drosophila_gpsim_fragment)

  # Get the target probe names
  library(annotate)
  aliasMapping <- getAnnMap("ALIAS2PROBE",
                    annotation(drosophila_gpsim_fragment))
  twi <- get('twi', env=aliasMapping)
  fbgnMapping <- getAnnMap("FLYBASE2PROBE",
                   annotation(drosophila_gpsim_fragment))
  targetProbe <- get('FBgn0035257', env=fbgnMapping)

  # Learn the model
  model <- GPLearn(drosophila_gpsim_fragment,
                   TF=twi, targets=targetProbe,
                   useGpdisim=TRUE, quiet=TRUE)

  # Plot it
  GPPlot(model, nameMapping=getAnnMap("FLYBASE",
                   annotation(drosophila_gpsim_fragment)))
}}
