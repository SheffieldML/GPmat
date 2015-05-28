function prior = normuniPriorParamInit(prior)

% NORMUNIPRIORPARAMINIT Normal uniform prior model's parameter initialisation.
% FORMAT
% DESC initialises the parameters of the normal uniform prior with some
% default parameters.
% ARG prior : prior structure to be initialised.
% RETURN prior : prior structure with initial values in place.
% 
% SEEALSO : priorCreate, gammaPriorParamInit, gaussianPriorParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% PRIOR

prior.sigma = 0.1;
prior.width = 2;
prior.transforms.index = [1 2];
prior.transforms.type = optimiDefaultConstraint('positive');
prior.nParams = 2;
