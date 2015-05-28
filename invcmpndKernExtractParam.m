function [params, names] = invcmpndKernExtractParam(kern)

% INVCMPNDKERNEXTRACTPARAM Extract parameters from the INVCMPND kernel structure.
% FORMAT
% DESC Extract parameters from the inv. precision compound kernel matrix into a vector of
% parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the
% kernel. The vector of 'transforms' is assumed to be empty
% here. Any transformations of parameters should be done in
% component kernels.
%
% SEEALSO cmpndKernExtractParam, invcmpndKernParamInit, invcmpndKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Andreas C. Damianou, 2012

% KERN


if nargout > 1
    [params, names] = cmpndKernExtractParam(kern);
else
    params = cmpndKernExtractParam(kern);
end
