function modelDisplay(model, varargin)

% MODELDISPLAY Display a text output of a model.

% MLTOOLS

% Check if the model has display code.
if exist([model.type 'Display'])==2
  fhandle = str2func([model.type 'Display']);
  fhandle(model, varargin{:});
end
