function dir = datasetsDirectory

% DATASETSDIRECTORY Returns directory where data is stored.

% DATASETS

% by default return the directory where this file is.
fullSpec = which('datasetsDirectory');
ind = max(find(fullSpec == filesep));
dir = fullSpec(1:ind);
