function kern = sdlfmvKernExpandParam(kern, params)

% SDLFMVKERNEXPANDPARAM Pass parameters from params to SDLFMV kernel
% FORMAT
% DESC returns a switching dynamical kernel structure for velocities filled 
% with the parameters in the given vector. This is used as a helper 
% function to enable parameters to be optimised in, for example, the NETLAB
% optimisation functions.
% ARG kern : the kernel structure in which the parameters are to be
% placed.
% ARG param : vector of parameters which are to be placed in the
% kernel structure.
% RETURN kern : kernel structure with the given parameters in the
% relevant locations.
%
% SEEALSO : sdlfmvKernParamInit, sdlfmvKernExtractParam, kernExpandParam
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

% This is just a wrapper function to call the expandParam function of the
% SDLFM structure.

kern = sdlfmKernExpandParam(kern, params);

