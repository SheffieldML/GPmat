function kern = rbfKernExpandParam(kern, params)

% RBFKERNEXPANDPARAM Create kernel structure from RBF kernel's parameters.
% FORMAT
% DESC returns a radial basis function kernel structure filled with the
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
% SEEALSO : rbfKernParamInit, rbfKernExtractParam, kernExpandParam
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% KERN

kern.inverseWidth = params(1);
kern.variance = params(2);
