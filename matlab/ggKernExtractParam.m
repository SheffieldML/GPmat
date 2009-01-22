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

% KERN

params = [kern.precision_u' kern.precision_y' kern.sigma2_u kern.sigma2_y kern.translation'];

if nargout > 1
    unames = cell(1,kern.inputDimension);
    ynames = cell(1,kern.inputDimension);
    mu_names = cell(1,kern.inputDimension);
    for i=1:kern.inputDimension,
        unames{i}=['inverse width latent (' num2str(i) ',' num2str(i) ')'];
        ynames{i}=['inverse width output (' num2str(i) ',' num2str(i) ')'];
        mu_names{i}=['mean (' num2str(i) ')'];
    end    
    names = {unames{:}, ynames{:}, 'variance latent', 'variance output' , mu_names{:}};
end