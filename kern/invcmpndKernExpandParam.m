function kern = invcmpndKernExpandParam(kern, params)


% INVCMPNDKERNEXPANDPARAM Create kernel structure from INVCMPND kernel's parameters.
% FORMAT
% DESC returns a inv-compound kernel structure filled with the
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
% SEEALSO : invcmpndKernParamInit, invcmpndKernExtractParam, cmpndkernExpandParam
%
% COPYRIGHT : Andreas C. Damianou, 2012

% KERN


kern = cmpndKernExpandParam(kern, params);
