function lvmResultsDynamic(modelType, dataSet, number, dataType, varargin)

% LVMRESULTSDYNAMIC Load a results file and visualise them.
% FORMAT
% DESC loads results of a latent variable model and visualises them.
% ARG modelType : the type of model ran on the data set.
% ARG dataSet : the name of the data set to load.
% ARG number : the number of the run used.
% ARG dataType : the type of data to visualise.
% ARG arg1, arg2, arg3 ... : additional arguments to be passed to the
% lvmVisualise command.
%
% SEEALSO : lvmLoadResult, lvmVisualise
%
% COPYRIGHT : Neil D. Lawrence, 2008
  
% MLTOOLS

[model, lbls] = lvmLoadResult(modelType, dataSet, number);

% Visualise the results
if size(model.X, 2) > 1
  lvmVisualise(model, lbls, [dataType 'Visualise'], [dataType 'Modify'], ...
                 varargin{:});
else  
  error('No visualisation code for data of this latent dimension.');
end
