function lvmResultsDynamic(modelType, dataSet, number, dataType, varargin)

% LVMRESULTSDYNAMIC Load a results file and visualise them.
%
%	Description:
%
%	LVMRESULTSDYNAMIC(MODELTYPE, DATASET, NUMBER, DATATYPE, ...) loads
%	results of a latent variable model and visualises them.
%	 Arguments:
%	  MODELTYPE - the type of model ran on the data set.
%	  DATASET - the name of the data set to load.
%	  NUMBER - the number of the run used.
%	  DATATYPE - the type of data to visualise.
%	  ... - additional arguments to be passed to the lvmVisualise
%	   command.
%	
%
%	See also
%	LVMLOADRESULT, LVMVISUALISE


%	Copyright (c) 2008 Neil D. Lawrence


[model, lbls] = lvmLoadResult(modelType, dataSet, number);

% Visualise the results
if size(model.X, 2) > 1
  lvmVisualise(model, lbls, [dataType 'Visualise'], [dataType 'Modify'], ...
                 varargin{:});
else  
  error('No visualisation code for data of this latent dimension.');
end