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
%
% MODIFICATIONS : Mauricio A. Alvarez, 2009.

% KERN

if kern.isArd
    params = kern.precisionT';
else
    params = kern.precisionT;
end
params(end+1) = kern.sigma2Noise;
if nargout > 1
    if kern.isArd
        unames = cell(size(kern.precisionT,1),1);
        for i=1:size(kern.precisionT,1),
            unames{i}=['VIK inverse width latent (' num2str(i) ',' num2str(i) ')'];
        end
        names = unames(:)';
        names = {names{:}, 'variance'};
    else
        unames = cell(numel(kern.precisionT),1);
        cont = 0;
        for j=1:size(kern.precisionT,2),
            cont = cont + 1;
            unames{cont}=['VIK ' num2str(j) ' inverse width latent'];
        end
        names = unames(:)';
        names = {names{:}, 'variance'};
    end
end