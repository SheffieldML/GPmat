function gpTimeDynamicsDisplay(model, varargin)

% GPTIMEDYNAMICSDISPLAY Display a GP time dynamics model.
% FORMAT
% DESC displays in human readable form the contents of the GP time dynamics
% model.
% ARG model : the model structure to be displaced.
% ARG spaceNum : number of spaces to place before displaying model
% structure.
%
% SEEALSO : gpDisplay, gpDynamicsCreate, modelDisplay.
%
% COPYRIGHT : Neil D. Lawrence, 2006

% FGPLVM

gpDisplay(model, varargin{:});
