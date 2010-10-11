\name{export.scores}
\Rdversion{1.0}
\alias{export.scores}
\title{Export results to an SQLite database}
\description{
  Exports the results to an SQLite database which can then be browsed
  with a result browser. The function will export log likelihoods,
  z-scores, model figures and gene aliases.
}
\usage{
  export.scores(scores, datasetName='', experimentSet='',
    databaseFile='database.sqlite', preprocData=NULL, models=NULL,
    figpath=NULL, aliasTypes=c("SYMBOL", "GENENAME", "ENTREZID"),
    datasetSource='', datasetDescription='',
    datasetSaveLocation='', datasetFigureFilename='',
    experimentTimestamp=as.character(Sys.Date()),
    figureDesc='', figurePrio=0, regulator=NULL)
}
\arguments{
  \item{scores}{The scoreList to export.}
  \item{datasetName}{Name of the dataset in the database.}
  \item{experimentSet}{Name of the experiment set in the database.}
  \item{databaseFile}{Filename of the database. New database is created if
    the file does not exist.}
  \item{preprocData}{Preprocessed data. This is required in order to generate
    models and figures and to calculate z-scores. Also, inserting aliases
    requires preprocessed data.}
  \item{models}{Learned models. If not given, the function will generate
    models if preprocessed data is available.}
  \item{figpath}{Figure path. If this is defined, the function will generate
    figures to the given path instead of inserting them to the database.}
  \item{aliasTypes}{Types of aliases that are inserted to the database.}
  \item{datasetSource}{Additional information that is inserted to the
    database if defined.}
  \item{datasetDescription}{Additional information that is inserted to the
    database if defined.}
  \item{datasetSaveLocation}{Additional information that is inserted to the
    database if defined.}
  \item{datasetFigureFilename}{Additional information that is inserted to the
    database if defined.}
  \item{experimentTimestamp}{Timestamp that is inserted to the database.
    The default value is current date in ISO-8601 format.}
  \item{figureDesc}{Additional information that is inserted to the
    database if defined.}
  \item{figurePrio}{Additional information that is inserted to the
    database if defined.}
  \item{regulator}{If defined, override the regulator name from scoreList.}
}
\author{Miika-Petteri Matikainen, Antti Honkela}
\seealso{
  \code{\link{GPRankTargets}, \link{GPRankTFs}}.
}
\examples{\dontrun{
  # Load a mmgmos preprocessed fragment of the Drosophila developmental
  # time series
  data(drosophila_gpsim_fragment)

  # FBgn names of target genes
  targets <- c('FBgn0003486', 'FBgn0033188', 'FBgn0035257')
  # Load gene annotations
  library(annotate)
  aliasMapping <- getAnnMap("ALIAS2PROBE",
                  annotation(drosophila_gpsim_fragment))

  # Get the probe identifier for TF 'twi'
  twi <- get('twi', env=aliasMapping)
  # Load alternative gene annotations
  fbgnMapping <- getAnnMap("FLYBASE2PROBE",
                 annotation(drosophila_gpsim_fragment))

  # Get the probe identifiers for target genes
  targetProbes <- mget(targets, env=fbgnMapping)

  # Rank the targets, filtering weakly expressed genes with average
  # expression z-score below 1.8
  scores <- GPRankTargets(drosophila_gpsim_fragment, TF=twi,
                          testTargets=targetProbes,
                          options=list(quiet=TRUE),
                          filterLimit=1.8)

  # Export data from scoreList and preprocessed data to a database
  export.scores(scores, datasetName='Drosophila',
                experimentSet='GPSIM/GPDISIM',
                database='database.sqlite',
                preprocData=drosophila_gpsim_fragment,
                aliasTypes=c('SYMBOL', 'GENENAME', 'FLYBASE', 'ENTREZID'))
}}
\keyword{export}
