function lfmResultsDynamic(dataSet, number, dataType, varargin)

% LFMRESULTSDYNAMIC Load a results file and visualise them.
% FORMAT
% DESC loads results of a latent variable model and visualises them.
% ARG dataSet : the name of the data set to load.
% ARG number : the number of the run used.
% ARG dataType : the type of data to visualise.
% ARG arg1, arg2, arg3 ... : additional arguments to be passed to the
% lvmVisualise command.
%
% SEEALSO : lfmVisualise
%
% COPYRIGHT : Neil D. Lawrence, 2008
%
% MODIFICATIONS : Mauricio Alvarez, 2009

% MLTOOLS

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
