function prior = gammaPriorParamInit(prior)

% GAMMAPRIORPARAMINIT Gamma prior model's parameter initialisation.
% FORMAT
% DESC initialises the parameters of the gamma prior with some
% default parameters.
% ARG prior : prior structure to be initialised.
% RETURN prior : prior structure with initial values in place.
% 
% SEEALSO : priorCreate
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% PRIOR

prior.a = 1e-6;
prior.b = 1e-6;

prior.transforms.index = [1 2];
prior.transforms.type = optimiDefaultConstraint('positive');
prior.nParams = 2;
