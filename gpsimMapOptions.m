function options = gpsimMapOptions(numGenes)

% GPSIMMAPOPTIONS Creates a set of default options for a GPSIMMAP model.
% FORMAT
% DESC returns a default options stucture for a GPSIMMAP model.
% ARG numGenes : the number of Genes that the model will consider.
% RETURN options : the options structure.
% 
% SEEALSO : gpsimMapCreate
%
% COPYRIGHT : Neil D. Lawrence, 2006

% SHEFFIELDML

options = gpsimOptions;

% Number of points for the numerical path integration (complexity
% is cubic in this value).
options.intPoints = 100;
options.kern = 'rbf';
options.nonLinearity = 'exp';

if nargin < 1
  options.B = [];
  options.D = [];
  options.S = [];
else
  options.B = ones(1, numGenes);
  options.D = ones(1, numGenes);
  options.S = ones(1, numGenes);
end
