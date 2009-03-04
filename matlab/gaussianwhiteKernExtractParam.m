function [params, names] = gaussianwhiteKernExtractParam(kern)

% GAUSSIANWHITEKERNEXTRACTPARAM Extract parameters from the gaussian white 
%                               kernel structure.
% FORMAT
% DESC extracts parameters from the
%	gaussian white kernel structure into a vector of parameters for
%	optimisation.
% RETURN param : vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel matrix, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
% ARG kern : the kernel structure containing the parameters to be
%	   extracted.
%
% DESC extracts parameters and
%	parameter names from the gaussian white kernel structure.
% RETURN param :vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel matrix, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
% RETURN names : cell array of strings giving names to the parameters.
% ARG kern : the kernel structure containing the parameters to be
%	   extracted.
%
% SEEALSO : gaussianwhiteKernParamInit, gaussianwhiteKernExpandParam,
% kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008

% KERN

params = kern.precisionT';
params(end+1) = kern.sigma2Noise;
if nargout > 1
    unames = cell(size(kern.precisionT,1),1);
    for i=1:size(kern.precisionT,1),
        unames{i}=['inverse width latent (' num2str(i) ',' num2str(i) ')'];
    end
    names = unames(:)';
    names = {names{:}, 'variance latent'};
end