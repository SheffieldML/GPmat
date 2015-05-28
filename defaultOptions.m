function options = defaultOptions;

% DEFAULTOPTIONS The default options for optimisation.
% FORMAT
% DESC returns a default options vector for optimisation.
% RETURN options : the default options vector.
% 
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
%
% SEEALSO : scg, conjgrad, quasinew

% NDLUTIL

options = [0,  1e-4, 1e-4, 1e-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1e-8, 0.1, 0];
