function [params, names] = lfmwhiteKernExtractParam(kern)

% LFMWHITEKERNEXTRACTPARAM Extract parameters from the LFM-WHITE kernel
% structure.
% FORMAT
% DESC extracts parameters from the LFM-White (Latent Force Model - White)
% kernel structure into a vector of parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel structure, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% FORMAT
% DESC extracts parameters and parameter names from the LFM-White (Latent
% Force Model - White) kernel structure.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If the
% field 'transforms' is not empty in the kernel structure, the parameters
% will be transformed before optimisation (for example positive only
% parameters could be logged before being returned).
% RETURN names : cell array of strings containing names for each
% parameter.
%
% SEEALSO lfmwhiteKernParamInit, lfmwhiteKernExpandParam, kernExtractParam,
% scg, conjgrad
%
% COPYRIGHT : David Luengo, 2009
%
% KERN


params = [kern.mass kern.spring kern.damper kern.variance kern.sensitivity];
if nargout > 1
  names = {'mass', 'spring', 'damper', 'variance', 'sensitivity'};
end
