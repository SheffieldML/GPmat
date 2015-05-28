function [params, names] = lfmaKernExtractParam(kern)

% LFMAKERNEXTRACTPARAM Extract parameters from the LFMA kernel structure.
% FORMAT
% DESC Extract parameters from the LFMA kernel structure
% into a vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. 
%
% FORMAT
% DESC Extract parameters from the LFMA kernel structure
% into a vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. 
% RETURN names : cell array of strings containing parameter names.
%
% SEEALSO lfmKernParamInit, lfmKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

[params, names] = lfmKernExtractParam(kern);
