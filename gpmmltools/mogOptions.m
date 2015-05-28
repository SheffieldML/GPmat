function options = mogOptions(numComp)

% MOGOPTIONS Sets the default options structure for MOG models.
% FORMAT
% DESC sets the default options structure for mixtures of
% Gaussians models.
% ARG numComponents : number of components in the mixture model.
% RETURN options : structure containing the default options.
%
% SEEALSO : mogCreate
% 
% COPYRIGHT : Neil D. Lawrence, 2006, 2008

% MLTOOLS

options.numComponents = numComp;
options.covtype = 'ppca';
% Whether it is an infinite mixture (false by default);
options.isInfinite = false;
options.a1=1;
