function lvmResultsClick(modelType, dataSet, number, dataType, varargin)

% LVMRESULTSCLICK Load a results file and visualise them with clicks
%
%	Description:
%
%	LVMRESULTSCLICK(MODELTYPE, DATASET, NUMBER, DATATYPE, ...) loads
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


%	Copyright (c) 2008, 2009 Neil D. Lawrence


[model, lbls] = lvmLoadResult(modelType, dataSet, number);

% Visualise the results
switch size(model.X, 2) 
 case 2
  lvmClickVisualise(model, lbls, [dataType 'Visualise'], [dataType 'Modify'], ...
                 varargin{:});
  
 otherwise 
  lvmClickVisualise(model, lbls, [dataType 'Visualise'], [dataType 'Modify'], ...
                 varargin{:});
end
