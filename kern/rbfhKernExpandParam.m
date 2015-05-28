function kern = rbfhKernExpandParam(kern, params)

% RBFHKERNEXPANDPARAM Create kernel structure from RBFH kernel's parameters.
% FORMAT
% DESC returns a radial basis function HEAT kernel structure filled with 
% the parameters in the given vector. This is used as a helper function to
% enable parameters to be optimised in, for example, the NETLAB
% optimisation functions.
% ARG kern : the kernel structure in which the parameters are to be
% placed.
% ARG param : vector of parameters which are to be placed in the
% kernel structure.
% RETURN kern : kernel structure with the given parameters in the
% relevant locations.
%
% SEEALSO : rbfhKernParamInit, rbfhKernExtractParam, kernExpandParam
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

kern.inverseWidthTime = params(1);
kern.inverseWidthSpace = params(2);
