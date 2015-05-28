function lfmResultsDynamic(dataSet, number, dataType, varargin)

% LFMRESULTSDYNAMIC Load a results file and visualise them.
%
%	Description:
%
%	LFMRESULTSDYNAMIC(DATASET, NUMBER, DATATYPE, ...) loads results of a
%	latent variable model and visualises them.
%	 Arguments:
%	  DATASET - the name of the data set to load.
%	  NUMBER - the number of the run used.
%	  DATATYPE - the type of data to visualise.
%	  ... - additional arguments to be passed to the lvmVisualise
%	   command.
%	
%	
%
%	See also
%	LFMVISUALISE


%	Copyright (c) 2008 Neil D. Lawrence


%	With modifications by Mauricio Alvarez 2009


%[model, lbls] = lvmLoadResult(modelType, dataSet, number);
capName = dataSet;
capName(1) = upper(capName(1));
load(['dem' capName num2str(number) '.mat'], 'model');
%load(['additional' capName '.mat'], 'skel', 'initPos');  

% Visualise the results
switch model.nlf
 case 2
  lfmVisualise(model, [dataType 'Visualise'], [dataType 'Modify'], varargin{:});
  
 otherwise 
  error('No visualisation code for data with this number of latent forces.');
end
