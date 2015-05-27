function [params, names] = lfmKernExtractParam(kern)

% LFMKERNEXTRACTPARAM Extract parameters from the LFM kernel structure.
% FORMAT
% DESC Extract parameters from the single input motif kernel structure
% into a vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. 
%
% FORMAT
% DESC Extract parameters and their names from the single input
% motif kernel structure.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. 
% RETURN names : cell array of strings containing parameter names.
%
% SEEALSO lfmKernParamInit, lfmKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2007

% KERN

params = [kern.mass, kern.spring, kern.damper,  kern.inverseWidth, kern.sensitivity];
if nargout > 1
  names = {'mass', 'spring', 'damper', 'inverse width', 'sensitivity'};
end
