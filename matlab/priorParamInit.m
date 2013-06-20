function prior = priorParamInit(prior, options)

% PRIORPARAMINIT Prior model's parameter initialisation.
% FORMAT
% DESC wrapper function for initialising prior distributions parameters.
% ARG prior : structure to initialise.
% RETURN prior : initialised prior structure.
%
% SEEALSO : priorCreate
%
% COPYRIGHT : Neil D. Lawrence, 2003
  
% PRIOR

if nargin > 1,
  prior = feval([prior.type 'PriorParamInit'], prior, options);
else
  prior = feval([prior.type 'PriorParamInit'], prior);
end
