function prior = mixturePriorParamInit(prior, options)

% MIXTUREPRIORPARAMINIT Mixture prior model's parameter initialisation.
% FORMAT
% DESC initialises the parameters of the mixture
% prior with some default parameters.
% ARG prior : prior structure to be initialised.
% ARG options : prior options, a structure with fields
% 'types': cell array of types; 'weights': vector of weights (optional)
% RETURN prior : prior structure with initial values in place.
% 
% SEEALSO : priorCreate
%
% COPYRIGHT : Antti Honkela, 2013

% SHEFFIELDML

if isfield(options, 'weights'),
  prior.weights = options.weights;
else
  prior.weights = 1/length(options.types) * ones(1, length(options.types));
end
prior.comp = cell(1, length(options.types));

prior.compNParams = zeros(size(prior.comp));
compIsBounded = zeros(size(prior.comp));

for k=1:length(options.types),
  prior.comp{k} = priorCreate(options.types{k});
  compNParams(k) = prior.comp{k}.nParams;
  compIsBounded(k) = isfield(prior.comp{k}, 'isBounded') && prior.comp{k}.isBounded;
end

prior.nParams = sum(compNParams);
if all(compIsBounded),
  prior.isBounded = 1;
elseif ~any(compIsBounded),
  prior.isBounded = 0;
else
  error('mixturePriorParamInit: mixing unbounded and bounded priors is not supported!');
end
