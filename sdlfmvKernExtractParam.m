function [params, names] = sdlfmvKernExtractParam(kern)

% SDLFMVKERNEXTRACTPARAM Extract parameters from the SDLFMV kernel structure.
% FORMAT
% DESC Extract parameters from the switching dynamical LFM kernel structure
% for the velocities into a vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. 
%
% FORMAT
% DESC Extract parameters and their names from the switching dynamical
% LFM kernel structure for the velocities.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. 
% RETURN names : cell array of strings containing parameter names.
%
% SEEALSO sdlfmKernExtractParam, sdlfmKernExpandParam, kernExtractParam 
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

% This is just a wrapper function that calls the extractParam function for
% the SDLFM kernel.

if nargout > 1
    [params, names] = sdlfmKernExtractParam(kern);
else
    params = sdlfmKernExtractParam(kern);
end
