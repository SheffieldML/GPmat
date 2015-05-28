function kern = dexpKernExpandParam(kern, params)

% DEXPKERNEXPANDPARAM Create a kernel structure from the double exponential
% kernel's parameters.
%
% FORMAT
% DESC returns a double exponential kernel kernel structure filled with the
% parameters in the given vector. This is used as a helper function to
% enable parameters to be optimised in, for example, the NETLAB
% optimisation functions.
% ARG kern : the kernel structure in which the parameters are to be
% placed.
% ARG params : vector of parameters which are to be placed in the
% kernel structure.
% RETURN kern : kernel structure with the given parameters in the
% relevant locations.
%
% SEEALSO : dexpKernParamInit, dexpKernExtractParam, kernExpandParam
%
% COPYRIGHT : David Luengo, 2009

% KERN


kern.decay = params(1);
kern.variance = params(2);
