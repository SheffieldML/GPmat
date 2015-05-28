function prior = laplacePriorParamInit(prior)

% LAPLACEPRIORPARAMINIT Laplace prior model's parameter initialisation.
% FORMAT
% DESC initialises the parameters of the Laplace prior with some
% default parameters.
% ARG prior : prior structure to be initialised.
% RETURN prior : prior structure with initial values in place.
% 
% SEEALSO : priorCreate, gammaPriorParamInit, gaussianPriorParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% PRIOR

prior.precision = 1;

prior.transforms.index = [1];
prior.transforms.type = optimiDefaultConstraint('positive');
prior.nParams = 1;
