function options = dnetOptions(mappingType, latentPoints, mappingOptions)

% DNETOPTIONS Options for a density network.
% FORMAT
% DESC returns the default options for a density network.
% RETURN options : default options structure for density network.
%
% FORMAT
% returns the default options for a density network given a
% number of hidden units.
% ARG mappingType : type of mapping to use in density network.
% ARG latentPoints : number of latent points to use. If give as a vector,
% layout is assumed to be a grid with rows and columns given by the vector.
% ARG mappingOptions : the options for the mapping model.
% RETURN options : default options structure for density network with the
% specified number of hidden units.
%
% SEEALSO : dnetCreate, modelOptions
%
% COPYRIGHT : Neil D. Lawrence, 2008

% MLTOOLS

if nargin < 1
  % This is the GTM settings
  mappingType = 'rbf';
end
if nargin < 2
  latentPoints = [16 16];
end
if nargin < 3
  options = modelOptions(mappingType);
end

options.mappingOptions = options;
options.mappingType = mappingType;
options.latentSamples = latentPoints;

options.alpha = 1;

options.activeFunc = 'linear';
options.optimiser = optimiDefaultOptimiser;

options.initX = 'ppca'; % Initialise X by ppca by default.
if length(latentPoints)>1
  options.M = prod(latentPoints);
  options.grid = latentPoints;
  options.basisStored = true;
else
  options.grid =[];
  options.M = latentPoints;
  options.basisStored = false;
end
