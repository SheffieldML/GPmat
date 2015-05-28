function prior = uniformPriorParamInit(prior)

% UNIFORMPRIORPARAMINIT Uniform prior model's parameter initialisation.
% FORMAT
% DESC initialises the parameters of the 2 parameter uniform
% prior with some default parameters.
% ARG prior : prior structure to be initialised.
% RETURN prior : prior structure with initial values in place.
% 
% SEEALSO : priorCreate
%
% COPYRIGHT : Antti Honkela, 2012

% PRIOR

prior.a = 0;
prior.b = 1;

prior.nParams = 0;
prior.isBounded = 1;
