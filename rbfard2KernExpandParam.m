function kern = rbfard2KernExpandParam(kern, params)

% RBFARD2KERNEXPANDPARAM Create kernel structure from RBFARD2 kernel's parameters.
% FORMAT
% DESC returns a automatic relevance determination radial basis
% function kernel structure filled with the parameters in the given
% vector. This is used as a helper function to enable parameters to
% be optimised in, for example, the NETLAB optimisation functions.
% ARG kern : the kernel structure in which the parameters are to be
% placed.
% ARG param : vector of parameters which are to be placed in the
% kernel structure.
% RETURN kern : kernel structure with the given parameters in the
% relevant locations.
%
% SEEALSO : rbfard2KernParamInit, rbfard2KernExtractParam, kernExpandParam
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
%
% COPYRIGHT : Michalis K. Titsias, 2009

% KERN

kern.variance = params(1);
kern.inputScales = params(2:end);
