function options = gpsimOptions

% GPSIMOPTIONS Creates a set of default options for a GPSIM model.
% FORMAT
% DESC returns a default options stucture for a GPSIM model.
% RETURN options : the options structure.
% 
% SEEALSO : gpsimCreate
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2008
%
% MODIFICATIONS : Pei Gao, 2008
  
% SHEFFIELDML

options.optimiser = 'conjgrad';
options.includeNoise = 0;
options.singleNoise = false;
