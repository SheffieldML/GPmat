function kern = ratquadKernExpandParam(kern, params)

% RATQUADKERNEXPANDPARAM Create kernel structure from RATQUAD kernel's parameters.
% FORMAT
% DESC returns a rational quadratic kernel structure filled with the
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
% SEEALSO : ratquadKernParamInit, ratquadKernExtractParam, kernExpandParam
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN

kern.alpha = params(1);
kern.lengthScale = params(2);
kern.variance = params(3);
