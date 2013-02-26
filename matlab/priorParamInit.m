function prior = priorParamInit(prior)

% PRIORPARAMINIT Prior model's parameter initialisation.
% FORMAT
% DESC wrapper function for initialising prior distributions parameters.
% ARG prior : structure to initialise.
% RETURN prior : initialised prior structure.
%
% SEEALSO : priorCreate
%
% COPYRIGHT : Neil D. Lawrence, 2003
  
% SHEFFIELDML

prior = feval([prior.type 'PriorParamInit'], prior);
