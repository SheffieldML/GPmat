function kern = expKernExpandParam(kern, params)

% EXPKERNEXPANDPARAM Create kernel structure from EXP kernel's parameters.
% FORMAT
% DESC returns a exponentiated kernel structure filled with the
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
% SEEALSO : expKernParamInit, expKernExtractParam, kernExpandParam
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN

kern.argument = kernExpandParam(kern.argument, params(2:end));
kern.variance = params(1);
% Deal with fact that white variance is exponentiated.
if isfield(kern.argument, 'whiteVariance')
  whiteVar = kern.argument.whiteVariance;
  kern.whiteVariance = kern.variance*(exp(2*whiteVar)-exp(whiteVar));
end
