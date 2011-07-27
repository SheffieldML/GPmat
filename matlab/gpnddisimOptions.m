function options = gpnddisimOptions

% GPNDDISIMOPTIONS Creates a set of default options for a GPNDDISIM model.
% FORMAT
% DESC returns a default options stucture for a GPNDDISIM model.
% RETURN options : the options structure.
% 
% SEEALSO : gpnddisimCreate
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2008
%
% COPYRIGHT : Pei Gao, 2008
%
% COPYRIGHT :  Jaakko Peltonen, 2011
  
% GPNDDISIM


options=struct();
options.includeNoise=1;
options.addPriors=0;
options.optimiser='quasinew';
options.use_disimstartmean=1;
options.use_fixedrnavariance=0;
options.fix=[];

