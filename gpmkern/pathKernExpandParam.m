function kern = pathKernExpandParam(kern, params)

% PATHKERNEXPANDPARAM Create kernel structure from PATH kernel's parameters.
% FORMAT
% DESC returns a path kernel structure filled with the
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
% SEEALSO : pathKernParamInit, pathKernExtractParam, kernExpandParam
%
% COPYRIGHT : Andrea Baisero, Carl Henrik Ek, 2013

% SHEFFIELDML


kern.cd=(params(1))^2;
kern.chv=params(2);
kern.gkern=kernExpandParam(kern.gkern,params(3:end));

kern = pathKernUpdateWMat(kern,[],true);