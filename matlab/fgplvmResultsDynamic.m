function fgplvmResultsDynamic(dataSet, number, dataType, varargin)

% FGPLVMRESULTSDYNAMIC Load a results file and visualise them.
%
% fgplvmResultsDynamic(dataSet, number, dataType, varargin)
%

% Copyright (c) 2006 Neil D. Lawrence
% fgplvmResultsDynamic.m version 1.1


  
[model, lbls] = fgplvmLoadResult(dataSet, number);

% Visualise the results
switch size(model.X, 2) 
 case 2
  fgplvmVisualise(model, lbls, [dataType 'Visualise'], [dataType 'Modify'], ...
                 varargin{:});
  
 otherwise 
  error('No visualisation code for data of this latent dimension.');
end