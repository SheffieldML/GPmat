function [params, names] = gaussianKernExtractParam(kern)

% GAUSSIANKERNEXTRACTPARAM Extract parameters from the gaussian kernel structure.
% FORMAT
% DESC extracts parameters from the
%	gaussian kernel structure into a vector of parameters for
%	optimisation.
% RETURN param : vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel matrix, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
% ARG kern : the kernel structure containing the parameters to be
%	   extracted.
%
% DESC extracts parameters and
%	parameter names from the gaussian kernel structure.
% RETURN param :vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel matrix, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
% RETURN names : cell array of strings giving names to the parameters.
% ARG kern : the kernel structure containing the parameters to be
%	   extracted.
%
% SEEALSO : gaussianKernParamInit, gaussianKernExpandParam,
% kernExtractParam, scg, conjgrad
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008
%
% MODIFICATIONS : Mauricio A. Alvarez, 2009

% KERN

params = kern.precisionU';
params(end+1) = kern.sigma2Latent;
if nargout > 1
    unames = cell(size(kern.precisionU,1),1);
    if exist([kern.type 'Names.txt'], 'file')
        fidNames = fopen([kern.type 'Names.txt'],'r');
        for i=1:size(kern.precisionU,1),           
            unames{i} = fgetl(fidNames);
        end
        fclose(fidNames);
    else
        for i=1:size(kern.precisionU,1),
            unames{i}=['inverse width latent ' num2str(i) '.'];
        end
        names = unames(:)';
        names = {names{:}, 'variance latent'};
    end
end
