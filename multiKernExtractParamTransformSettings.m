function [paramtransformsettings, names] = multiKernExtractParamTransformSettings(kern)

% MULTIKERNEXTRACTPARAMTRANSFORMSETTINGS Extract parameter transform settings 
% from the MULTI kernel structure.
% FORMAT
% DESC Extract parameters from the multiple output block kernel
% structure into a vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the
% kernel. The vector of transforms is assumed to be empty here, any
% transormation of parameters is assumed to be done in the
% component kernels.
%
% FORMAT
% DESC Extract parameters and parameter names from the multiple
% output block kernel structure.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the
% kernel. The vector of transforms is assumed to be empty here, any
% transormation of parameters is assumed to be done in the
% component kernels.
% RETURN names : cell array of strings containing parameter names.
%
% SEEALSO multiKernParamInit, multiKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN



if nargout > 1
  [paramtransformsettings, names] = cmpndKernExtractParamTransformSettings(kern);
else
  paramtransformsettings = cmpndKernExtractParamTransformSettings(kern);
end
