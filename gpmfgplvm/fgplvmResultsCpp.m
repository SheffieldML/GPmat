function fgplvmResultsCpp(fileName, dataType, varargin)

% FGPLVMRESULTSCPP Load a results file and visualise them.
% FORMAT
% DESC loads a model stored in the file name and visualizes results with
% it.
% ARG fileName : file name of the stored model.
% ARG dataType : type of visualization to perform.
% ARG opt1, opt2, ... : additional options to pass to the visualization.
%
% SEEALSO gplvmResultsCpp, lvmVisualise
%
% COPYRIGHT : Neil D. Lawrence, 2009
  
% FGPLVM
  
[model, lbls] = fgplvmReadFromFile(fileName);

lvmVisualise(model, lbls, [dataType 'Visualise'], [dataType 'Modify'], ...
             varargin{:});
