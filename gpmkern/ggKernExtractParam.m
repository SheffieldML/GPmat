function [params, names] = ggKernExtractParam(kern)

% GGKERNEXTRACTPARAM Extract parameters from the GG kernel structure.
% FORMAT
% DESC Extract parameters from the gaussian
%	gaussian kernel structure into a vector of parameters for optimisation.
% RETURN param : vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel matrix, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
% ARG kern : the kernel structure containing the parameters to be
%	   extracted.
%
% FORMAT
% DESC extract parameters and
%	their names from the gaussian gaussian kernel structure.
% RETURN param : vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel matrix, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
% RETURN names : cell array of strings containing parameter names.
% ARG kern : the kernel structure containing the parameters to be
%	   extracted.
%	
% SEEALSO : ggKernParamInit, ggKernExpandParam, kernExtractParam, scg,
% conjgrad
% 
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008
%
% MODIFICATIONS : Mauricio A. Alvarez, 2009

% KERN

params = [kern.precisionU(:)' kern.precisionG(:)' ...
    kern.sigma2Latent kern.sensitivity];

if nargout > 1    
    unames = cell(1,size(kern.precisionU,1));
    ynames = cell(1,size(kern.precisionG,1));
    if exist([kern.type 'Names.txt'], 'file')
        fidNames = fopen([kern.type 'Names.txt'],'r');
        for i=1:size(kern.precisionU,1),
            unames{i} = fgetl(fidNames);
            ynames{i} = fgetl(fidNames);
        end
        fclose(fidNames);
    else
        for i=1:size(kern.precisionU,1),
            unames{i}=['inverse width latent ' num2str(i) '.'];
            ynames{i}=['inverse width output ' num2str(i) '.'];
        end
    end
    names = {unames{:}, ynames{:}, 'variance latent', 'sensitivity'};
end
