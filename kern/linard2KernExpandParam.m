function kern = linard2KernExpandParam(kern, params)


% LINARD2KERNEXPANDPARAM Create kernel structure from LINARD2 kernel's parameters.
% FORMAT
% DESC returns a automatic relevance determination linear kernel structure filled with the
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
% SEEALSO : linard2KernParamInit, linard2KernExtractParam, kernExpandParam
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
%
% COPYRIGHT : Michalis K. Titsias, 2009

% KERN

kern.inputScales = params(1:end);
