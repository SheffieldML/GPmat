function kern = ouKernExpandParam(kern, params)

% OUKERNEXPANDPARAM Create kernel structure from OU kernel's parameters
% (see ouKernCompute or ouKernParamInit for a more detailed description of
% the OU kernel).
% FORMAT
% DESC returns a Ornstein-Uhlenbeck kernel kernel structure filled with the
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
% SEEALSO : ouKernParamInit, ouKernExtractParam, kernExpandParam
%
% COPYRIGHT : David Luengo, 2009

% KERN


kern.decay = params(1);
kern.variance = params(2);
