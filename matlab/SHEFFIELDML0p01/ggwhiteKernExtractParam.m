function [params, names] = ggwhiteKernExtractParam(kern)

% GGWHITEKERNEXTRACTPARAM Extract parameters from the GG WHITE kernel structure.
%
%	Description:
%
%	PARAM = GGWHITEKERNEXTRACTPARAM(KERN) Extract parameters from the
%	gaussian gaussian white kernel structure into a vector of parameters
%	for optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel matrix, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = GGWHITEKERNEXTRACTPARAM(KERN) extract parameters
%	and their names from the gaussian gaussian white kernel structure.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel matrix, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	  NAMES - cell array of strings containing parameter names.
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%	conjgrad
%	
%	
%
%	See also
%	GGWHITEKERNPARAMINIT, GGWHITEKERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, 


%	Copyright (c) 2008 Mauricio A. Alvarez and Neil D. Lawrence


%	With modifications by Mauricio A Alvarez 2009


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