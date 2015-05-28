function kern = gibbsKernExpandParam(kern, params)

% GIBBSKERNEXPANDPARAM Create kernel structure from GIBBS kernel's parameters.
% FORMAT
% DESC returns a Mark Gibbs's non-stationary kernel structure filled with the
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
% SEEALSO : gibbsKernParamInit, gibbsKernExtractParam, kernExpandParam
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN

kern.lengthScaleFunc = modelExpandParam(kern.lengthScaleFunc, params(1:end-1));
kern.variance = params(end);
