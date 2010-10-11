.openConnection <- function(dbPath) {
  drv <- SQLite()
  db <- dbConnect(drv, dbPath)
  return(db)
}

.createTables <- function(db) {
  tables <- cbind('CREATE TABLE genes (gene_id INTEGER PRIMARY KEY, probe_name VARCHAR,
CONSTRAINT genes_u UNIQUE (probe_name) ON CONFLICT IGNORE)',
   'CREATE TABLE gene_aliases (gene_id INTEGER, alias_id INTEGER, alias VARCHAR,
CONSTRAINT gene_aliases_u UNIQUE (gene_id, alias_id, alias) ON CONFLICT IGNORE,
CONSTRAINT gene_aliases_g_f FOREIGN KEY (gene_id) REFERENCES genes(gene_id),
CONSTRAINT gene_aliases_a_f FOREIGN KEY (alias_id) REFERENCES alias_annotation(alias_id))',
   'CREATE TABLE z_scores (gene_id INTEGER, dataset_id INTEGER, mean_z_score DOUBLE,
CONSTRAINT z_scores_u UNIQUE (gene_id, dataset_id) ON CONFLICT IGNORE,
CONSTRAINT z_scores_g_f FOREIGN KEY (gene_id) REFERENCES genes(gene_id),
CONSTRAINT z_scores_g_d FOREIGN KEY (dataset_id) REFERENCES expression_dataset_annotation(dataset_id))',
   'CREATE TABLE results (gene_id INTEGER, experiment_id INTEGER, log_likelihood DOUBLE, baseline_log_likelihood DOUBLE, params BLOB,
CONSTRAINT results_u UNIQUE (gene_id, experiment_id) ON CONFLICT REPLACE,
CONSTRAINT results_g_f FOREIGN KEY (gene_id) REFERENCES genes(gene_id),
CONSTRAINT results_e_f FOREIGN KEY (experiment_id) REFERENCES experiment_annotation(experiment_id))',
   'CREATE TABLE supplementary_data (gene_id INTEGER, supp_dataset_id INTEGER, value DOUBLE,
CONSTRAINT supplementary_data_u UNIQUE (gene_id, supp_dataset_id) ON CONFLICT IGNORE,
CONSTRAINT supplementary_data_g_f FOREIGN KEY (gene_id) REFERENCES genes(gene_id),
CONSTRAINT supplementary_data_s_f FOREIGN KEY (supp_dataset_id) REFERENCES supplementary_dataset_annotation(supp_dataset_id))',
   'CREATE TABLE target_sets (experiment_id INTEGER, gene_id INTEGER,
CONSTRAINT target_sets_u UNIQUE (experiment_id, gene_id) ON CONFLICT IGNORE,
CONSTRAINT target_sets_e_f FOREIGN KEY (experiment_id) REFERENCES experiment_annotation(experiment_id),
CONSTRAINT target_sets_g_f FOREIGN KEY (gene_id) REFERENCES genes(gene_id))',
   'CREATE TABLE figures (experiment_id INTEGER, figure_id INTEGER,
CONSTRAINT figures_u UNIQUE (experiment_id, figure_id) ON CONFLICT IGNORE,
CONSTRAINT figures_e_f FOREIGN KEY (experiment_id) REFERENCES experiment_annotation(experiment_id),
CONSTRAINT figures_f_f FOREIGN KEY (figure_id) REFERENCES figure_annotation(figure_id))',
   'CREATE TABLE figuredata (figure_id INTEGER, gene_id INTEGER, data BLOB,
CONSTRAINT figuredata_u UNIQUE (figure_id, gene_id) ON CONFLICT REPLACE,
CONSTRAINT figuredata_f_f FOREIGN KEY (figure_id) REFERENCES figure_annotation(figure_id),
CONSTRAINT figures_g_f FOREIGN KEY (gene_id) REFERENCES genes(gene_id))',
   'CREATE TABLE regulators (regulator_id INTEGER PRIMARY KEY, gene_id INTEGER, dataset_id INTEGER, regulator_name VARCHAR,
CONSTRAINT regulators_u UNIQUE (dataset_id, regulator_name) ON CONFLICT IGNORE,
CONSTRAINT regulators_g_f FOREIGN KEY (gene_id) REFERENCES genes(gene_id),
CONSTRAINT regulators_d_f FOREIGN KEY (dataset_id) REFERENCES expression_dataset_annotation(dataset_id))',
   'CREATE TABLE figure_annotation (figure_id INTEGER PRIMARY KEY, filename, name VARCHAR, description VARCHAR, priority INTEGER,
CONSTRAINT figure_annotation_u UNIQUE (filename, name) ON CONFLICT IGNORE)',
   'CREATE TABLE experiment_annotation (experiment_id INTEGER PRIMARY KEY, dataset_id INTEGER, regulator_id INTEGER, loop_variable INTEGER, model_translation BOOLEAN, number_of_parameters INTEGER, parameter_names VARCHAR, producer VARCHAR, timestamp VARCHAR, description VARCHAR, name VARCHAR,
CONSTRAINT experiment_annotation_u UNIQUE (dataset_id, regulator_id, loop_variable, model_translation, name) ON CONFLICT IGNORE,
CONSTRAINT experiment_annotation_d_f FOREIGN KEY (dataset_id) REFERENCES expression_dataset_annotation(dataset_id),
CONSTRAINT experiment_annotation_r_f FOREIGN KEY (regulator_id) REFERENCES regulators(regulator_id))',
   'CREATE TABLE experiment_set (experiment_set_id INTEGER PRIMARY KEY, parent_id INTEGER, name VARCHAR,
CONSTRAINT experiment_set_u UNIQUE (name) ON CONFLICT IGNORE,
CONSTRAINT experiment_set_p_f FOREIGN KEY (parent_id) REFERENCES experiment_set(experiment_set_id))',
   'CREATE TABLE experiment_set_experiments (experiment_set_id INTEGER, experiment_id INTEGER,
CONSTRAINT experiment_set_experiments_u UNIQUE (experiment_set_id, experiment_id) ON CONFLICT IGNORE,
CONSTRAINT experiment_set_experiments_s_f FOREIGN KEY (experiment_set_id) REFERENCES experiment_set(experiment_set_id),
CONSTRAINT experiment_set_experiments_e_f FOREIGN KEY (experiment_id) REFERENCES experiment_annotation(experiment_id))',
   'CREATE TABLE expression_dataset_annotation (dataset_id INTEGER PRIMARY KEY, dataset_name VARCHAR, species VARCHAR, source VARCHAR, platform VARCHAR, description VARCHAR, save_location VARCHAR, figure_filename VARCHAR,
CONSTRAINT expression_dataset_annotation_u UNIQUE (dataset_name) ON CONFLICT IGNORE)',
   'CREATE TABLE alias_annotation (alias_id INTEGER PRIMARY KEY, dataset_id INTEGER, alias_class VARCHAR, source VARCHAR, description VARCHAR,
CONSTRAINT alias_annotation_u UNIQUE (dataset_id, alias_class) ON CONFLICT IGNORE,
CONSTRAINT alias_annotation_d_f FOREIGN KEY (dataset_id) REFERENCES expression_dataset_annotation(dataset_id))',
   'CREATE TABLE supplementary_dataset_annotation (supp_dataset_id INTEGER PRIMARY KEY, supp_dataset_type INTEGER, supp_dataset_name VARCHAR, regulator_id INTEGER, source VARCHAR, platform VARCHAR, description VARCHAR,
CONSTRAINT supplementary_dataset_annotation_u UNIQUE (supp_dataset_name, regulator_id) ON CONFLICT IGNORE,
CONSTRAINT supplementary_dataset_annotation_r_f FOREIGN KEY (regulator_id) REFERENCES regulators(regulator_id))')

  VERSION <- 1

  tableCount <- length(dbListTables(db))
  if (tableCount > 0) {
    return(FALSE)
  }

  for (i in seq_along(tables)) {
    table <- tables[i]
    dbGetQuery(db, table)
  }

  dbGetQuery(db, paste("PRAGMA user_version = ", VERSION))
  return(TRUE)
}

.getVersion <- function(db) {
  return(dbGetQuery(db, "PRAGMA user_version")[1,1])
}

.addAndGetAliasId <- function(db, datasetId, aliasClass, source='', description='') {
  query <- paste("INSERT INTO alias_annotation",
                 "VALUES (null, @did, @class, @source, @desc)")
  dbGetPreparedQuery(db, query, data.frame(did = datasetId, class = aliasClass, source = source, desc = description))

  alias_query <- paste("SELECT aa.alias_id",
                       "FROM alias_annotation AS aa",
                       "WHERE aa.alias_class = @ac AND aa.dataset_id = @did")
  alias_id <- dbGetPreparedQuery(db, alias_query, data.frame(ac = aliasClass, did = datasetId))
  return(alias_id[1,1])
}

.addAndGetProbeGeneIds <- function(db, probeNames) {
  insert_query <- paste("INSERT INTO genes",
                        "VALUES (null, @probe)")

  data <- as.data.frame(probeNames)
  colnames(data) <- "probe"

  dbBeginTransaction(db)
  dbGetPreparedQuery(db, insert_query, data)
  dbCommit(db)

  query <- paste("SELECT g.probe_name, g.gene_id",
                 "FROM genes AS g",
                 "WHERE g.probe_name IN (@probe_names)")
  result <- dbGetPreparedQuery(db, query, data.frame(probe_names = probeNames))
  return(result)
}

.getProbeGeneId <- function(db, probe) {
  query <- paste("SELECT g.gene_id",
                 "FROM genes AS g",
                 "WHERE g.probe_name = @probe")
  return(dbGetPreparedQuery(db, query, data.frame(probe = probe)))
}

.addAndGetRegulatorId <- function(db, regulatorName, datasetId) {
  insert_query <- paste("INSERT INTO regulators",
                        "VALUES (null, @gene_id, @dataset_id, @regulator_name)")

  gene_id <- .getProbeGeneId(db, regulatorName)
  if (nrow(gene_id) == 0) {
    gene_id <- NA
  } else {
    gene_id <- gene_id[1,1]
  }

  dbGetPreparedQuery(db, insert_query, data.frame(gene_id = gene_id, dataset_id = datasetId, regulator_name = regulatorName))

  query <- paste("SELECT r.regulator_id",
                 "FROM regulators AS r",
                 "WHERE r.regulator_name = @regulator_name")
  return(dbGetPreparedQuery(db, query, data.frame(regulator_name = regulatorName)))
}

.addAndGetDatasetId <- function(db, name, species='', source='', platform='', description='', saveLocation='', figureFilename='') {
  query <- paste("INSERT INTO expression_dataset_annotation",
                 "VALUES (null, @name, @species, @source, @platform, @desc, @save_location, @figure_filename)")
  dbGetPreparedQuery(db, query, data.frame(name = name, species = species, source = source, platform = platform, desc = description, save_location = saveLocation, figure_filename = figureFilename))

  dataset_query <- paste("SELECT d.dataset_id",
                         "FROM expression_dataset_annotation AS d",
                         "WHERE d.dataset_name = @dname")
  dataset_id <- dbGetPreparedQuery(db, dataset_query, data.frame(dname = name))
  return(dataset_id[1,1])
}

.addAndGetExperimentId <- function(db, name, description, datasetId, regulatorId=NA, loopVariable=2, modelTranslation=FALSE, numberOfParameters=NA, parameterNames=NA, producer='', timestamp='') {
  insert_query <- paste("INSERT INTO experiment_annotation",
                        "VALUES (null, @dataset_id, @regulator_id, @loop_variable, @model_translation, @number_of_parameters, @parameter_names, @producer, @timestamp, @description, @name)")
  dbGetPreparedQuery(db, insert_query, data.frame(dataset_id = datasetId, regulator_id = regulatorId, loop_variable = loopVariable, model_translation = modelTranslation, number_of_parameters = numberOfParameters, parameter_names = parameterNames, producer = producer, timestamp = timestamp, description = description, name = name))

  query <- paste("SELECT e.experiment_id",
                 "FROM experiment_annotation AS e",
                 "WHERE e.name = @name AND e.dataset_id = @dataset_id AND e.regulator_id = @regulator_id AND e.loop_variable = @loop_variable AND e.model_translation = @model_translation")
  experiment_id <- dbGetPreparedQuery(db, query, data.frame(name = name, dataset_id = datasetId, regulator_id = regulatorId, loop_variable = loopVariable, model_translation = modelTranslation))
  return(experiment_id[1,1])
}

.addAndGetFigureId <- function(db, filename, name='', description='', priority=0) {
  insert_query <- paste("INSERT INTO figure_annotation",
                        "VALUES (null, @filename, @name, @description, @priority)")
  dbGetPreparedQuery(db, insert_query, data.frame(filename = filename, name = name, description = description, priority = priority))

  query <- paste("SELECT f.figure_id",
                 "FROM figure_annotation AS f",
                 "WHERE f.filename = @filename AND f.name = @name")
  figure_id <- dbGetPreparedQuery(db, query, data.frame(filename = filename, name = name))
  return(figure_id[1,1])
}

.addAndGetSupplementaryDataId <- function(db, name, regulatorId=NA, type=0, source='', platform='', description='') {
  insert_query <- paste("INSERT INTO supplementary_dataset_annotation",
                        "VALUES (null, @type, @name, @regulator_id, @source, @platform, @description)")
  dbGetPreparedQuery(db, insert_query, data.frame(type = type, name = name, regulator_id = regulatorId, source = source, platform = platform, description = description))

  query <- paste("SELECT s.supp_dataset_id",
                 "FROM supplementary_dataset_annotation AS s",
                 "WHERE s.name = @name AND s.regulator_id = @regulator_id")
  supp_id <- dbGetPreparedQuery(db, query, data.frame(name = name, regulator_id = regulatorId))
  return(supp_id[1,1])
}

.addAndGetExperimentSetId <- function(db, name, parentId) {
  insert_query <- paste("INSERT INTO experiment_set",
                        "VALUES (null, @parent_id, @name)")
  dbGetPreparedQuery(db, insert_query, data.frame(parent_id = parentId, name = name))

  query <- paste("SELECT s.experiment_set_id",
                 "FROM experiment_set AS s",
                 "WHERE s.name = @name")
  set_id <- dbGetPreparedQuery(db, query, data.frame(name = name, parent_id = parentId))
  return(set_id[1,1])
}

.addExperimentSetExperiments <- function(db, setId, experimentIds) {
  insert_query <- paste("INSERT INTO experiment_set_experiments",
                        "VALUES (@set_id, @experiment_id)")
  data <- as.data.frame(cbind(setId, experimentIds))
  names(data) <- cbind("set_id", "experiment_id")
  dbGetPreparedQuery(db, insert_query, data)
}

.addSupplementaryData <- function(db, suppDatasetId, suppData) {
  insert_query <- paste("INSERT INTO supplementary_data",
                        "VALUES (@gene_id, @supp_id, @value)")

  probes <- names(suppData)
  probe_gene_ids <- .addAndGetProbeGeneIds(db, probes)
  data <- as.data.frame(cbind(probe_gene_ids, suppDatasetId, suppData))
  data <- data[2:4] # remove probe_name column
  names(data) <- cbind("gene_id", "supp_id", "value")
  dbGetPreparedQuery(db, insert_query, data)
}

.addFigures <- function(db, experimentId, filename, name='', description='', priority=0, figureData=NULL) {
  insert_query <- paste("INSERT INTO figures",
                        "VALUES (@experiment_id, @figure_id)")

  figure_id <- .addAndGetFigureId(db, filename, name, description, priority)
  dbGetPreparedQuery(db, insert_query, data.frame(experiment_id = experimentId, figure_id = figure_id))

  if (length(figureData) > 0) {
    insert_query <- paste("INSERT INTO figuredata",
                          "VALUES (@figure_id, @gene_id, @data)")

    probes <- names(figureData)
    probe_gene_ids <- .addAndGetProbeGeneIds(db, probes)
    data <- as.data.frame(cbind(figure_id, probe_gene_ids, I(figureData)))
    data <- data[c(1,3,4)] # remove probe_name column

    # replace filename column by raw file contents
    for (k in seq_along(data[[3]])) {
      name <- data[[3]][[k]]
      to.read <- file(name, "rb")
      fileSize <- file.info(name)$size
      data[[3]][[k]] <- I(readBin(to.read, "raw", n=fileSize))
      close(to.read)
    }

    names(data) <- cbind("figure_id", "gene_id", "data")
    dbGetPreparedQuery(db, insert_query, data)
  }
}

.addResults <- function(db, experimentId, logLikelihoods, baselineLogLikelihoods=NA, params=NA) {
  insert_query <- paste("INSERT INTO results",
                        "VALUES (@gene_id, @experiment_id, @log_likelihood, @baseline_log_likelihood, @params)")
  probes <- names(logLikelihoods)
  probe_gene_ids <- .addAndGetProbeGeneIds(db, probes)
  data <- as.data.frame(cbind(probe_gene_ids, experimentId, logLikelihoods, baselineLogLikelihoods, I(lapply(params, function(x) serialize(x, NULL, ascii = FALSE)))))
  data <- data[2:6] # remove probe_name column
  names(data) <- cbind("gene_id", "experiment_id", "log_likelihood", "baseline_log_likelihood", "params")
  dbGetPreparedQuery(db, insert_query, data)
}

.addGeneAliases <- function(db, aliasId, geneId, aliasList) {
  insert_query <- paste("INSERT INTO gene_aliases",
                        "VALUES (@gene_id, @alias_id, @alias)")
  for (alias in aliasList) {
    dbGetPreparedQuery(db, insert_query, data.frame(gene_id = geneId, alias_id = aliasId, alias = alias))
  }
}

.addZscores <- function(db, datasetId, zscores) {
  insert_query <- paste("INSERT INTO z_scores",
                        "VALUES (@gene_id, @dataset_id, @zscore)")
  probes <- names(zscores)
  probe_gene_ids <- .addAndGetProbeGeneIds(db, probes)
  data <- as.data.frame(cbind(probe_gene_ids, datasetId, zscores))
  data <- data[2:4] # remove probe_name column
  names(data) <- cbind("gene_id", "dataset_id", "zscore")
  dbGetPreparedQuery(db, insert_query, data)
}

.addTargetSet <- function(db, experimentId, probes) {
  insert_query <- paste("INSERT INTO target_sets",
                        "VALUES (@experiment_id, @gene_id)")
  probe_gene_ids <- .addAndGetProbeGeneIds(db, probes)
  data <- as.data.frame(cbind(experimentId, probe_gene_ids))
  data <- data[c(1,3)] # remove probe_name column
  names(data) <- cbind("experiment_id", "gene_id")
  dbGetPreparedQuery(db, insert_query, data)
}

.addAliases <- function(db, probes, preprocData, datasetId, aliasType) {
  aliasSource <- annotation(preprocData)
  aliasDescription <- aliasType
  aliasClass <- sub("2PROBE", "", aliasType, ignore.case = TRUE)
  alias_id <- .addAndGetAliasId(db, datasetId, aliasClass, aliasSource, aliasDescription)
  aliasMapping <- getAnnMap(aliasType, annotation(preprocData))

  if (length(grep("2PROBE", aliasType, ignore.case = TRUE)) > 0) { # reversed mapping
    aliases <- revmap(aliasMapping)
  } else {
    aliases <- aliasMapping
  }

  aliases <- as.list(aliases)
  probe_gene_ids <- .addAndGetProbeGeneIds(db, probes)

  dbBeginTransaction(db)
  for (i in seq_along(probe_gene_ids[[1]])) {
    probe <- probe_gene_ids[i, 1]
    gene_id <- probe_gene_ids[i, 2]
    if (probe %in% names(aliases)) {
      alias_list <- aliases[probe]
      .addGeneAliases(db, alias_id, gene_id, alias_list)
    }
  }
  dbCommit(db)
}

.getLoopVariables <- function(db) {
  query <- paste("SELECT e.loop_variable",
                 "FROM experiment_annotation AS e")

  loopVariables <- dbGetQuery(db, query)$loop_variable
  return(loopVariables)
}


