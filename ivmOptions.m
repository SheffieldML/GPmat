function options = ivmOptions(varargin)

% IVMOPTIONS Return default options for IVM model.
% FORMAT
% DESC returns the default options in a structure for a IVM model.
% RETURN options : structure containing the default options for the
% given approximation type.
%
% SEEALSO : ivmCreate
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2005

% IVM

% bog-standard kernel.
options.kern = {'rbf', 'bias', 'white'};
options.numActive = 100;
options.noise = 'probit';
options.selectionCriterion = 'entropy';
options.display = 0;
options.kernIters = 100;
options.noiseIters = 0;
options.extIters = 8;

