function kern = polyKernExpandParam(kern, params)


% POLYKERNEXPANDPARAM Create kernel structure from POLY kernel's parameters.
% FORMAT
% DESC returns a polynomial kernel structure filled with the
% parameters in the given vector. This is used as a helper function to
% enable parameters to be optimised in, for example, the NETLAB
% optimisation functions.
% ARG kern : the kernel structure in which the parameters are to be
% placed.
% ARG param : vector of parameters which are to be placed in the
% kernel structure.
% RETURN kern : kernel structure with the given parameters in the
% relevant locations.
%
% SEEALSO : polyKernParamInit, polyKernExtractParam, kernExpandParam
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% KERN


%/~
if any(isinf(params))
  warning('params are infinite')
end
%~/

kern.weightVariance = params(1);
kern.biasVariance = params(2);
kern.variance = params(3);
