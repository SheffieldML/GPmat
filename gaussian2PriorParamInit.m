function prior = gaussian2PriorParamInit(prior, params)

% GAUSSIAN2PRIORPARAMINIT Gaussian prior model's parameter initialisation.
% This is the same as the gaussian prior but here the mean is not
% necessarily zero. Also, initial params are given optionally.
% FORMAT
% DESC initialises the parameters of the Gaussian prior with some
% default parameters.
% ARG prior : prior structure to be initialised.
% RETURN prior : prior structure with initial values in place.
% 
% SEEALSO : priorCreate, gammaPriorParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
% MODICICATIONS: Andreas Damianou, Carl Henrik Ek, 2013

% SHEFFIELDML

if nargin < 2
    precision = 1; mean = 0;
else
    precision = params(1);
    if length(params) > 1
        mean = params(2:end);
    else
        mean = 0;
    end
end

prior.precision = precision; % Assume isotropic variance if dim > 1
prior.mean = mean;
prior.dim = length(mean); 
prior.nParams = prior.dim + 1;
prior.transforms.index = 1:prior.nParams;
prior.transforms.type = optimiDefaultConstraint('positive');
