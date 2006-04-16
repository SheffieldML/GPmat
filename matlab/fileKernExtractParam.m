function params = fileKernExtractParam(kern)


% FILEKERNEXTRACTPARAM Extract parameters from the FILE kernel structure.
% FORMAT
% DESC Extract parameters from the stored file kernel matrix into a vector of
% parameters for optimisation.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the kernel. If
% the field 'transforms' is not empty in the kernel matrix, the
% parameters will be transformed before optimisation (for example
% positive only parameters could be logged before being returned).
%
% SEEALSO fileKernParamInit, fileKernExpandParam, kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
%
% KERN


params = [kern.variance];
