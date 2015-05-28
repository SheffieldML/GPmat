function [params, names] = ggwhiteKernExtractParam(kern)

% GGWHITEKERNEXTRACTPARAM Extract parameters from the GG WHITE kernel structure.
% FORMAT
% DESC Extract parameters from the gaussian
%	gaussian white kernel structure into a vector of parameters for optimisation.
% RETURN param : vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel matrix, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
% ARG kern : the kernel structure containing the parameters to be
%	   extracted.
%
% FORMAT
% DESC extract parameters and
%	their names from the gaussian gaussian white kernel structure.
% RETURN param : vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel matrix, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
% RETURN names : cell array of strings containing parameter names.
% ARG kern : the kernel structure containing the parameters to be
%	   extracted.
%	
% SEEALSO : ggwhiteKernParamInit, ggwhiteKernExpandParam, kernExtractParam, scg,
% conjgrad
% 
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008
%
% MODIFICATIONS : Mauricio A Alvarez, 2009

% KERN

params = [kern.precisionG' kern.sigma2Noise kern.variance];

if nargout > 1
    if kern.isArd
        ynames = cell(1,kern.inputDimension);
        for i=1:kern.inputDimension,
            ynames{i}=['inverse width output ' num2str(i) '.'];
        end
        names = {ynames{:}, 'variance latent', 'sensitivity'};
    else
        names = {'inverse width output 1.' , 'variance latent', 'sensitivity'};
    end    
end
