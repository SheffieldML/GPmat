function [params, names] = heatKernExtractParam(kern)

% HEATKERNEXTRACTPARAM Extract parameters from the HEAT kernel structure.
%
%	Description:
%
%	PARAM = HEATKERNEXTRACTPARAM(KERN) Extract parameters from the heat
%	kernel structure into a vector of parameters for optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel matrix, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = HEATKERNEXTRACTPARAM(KERN) Extract parameters and
%	their names from the heat kernel structure.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel matrix, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	  NAMES - cell array of strings containing parameter names.
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%	
%	
%
%	See also
%	% SEEALSO HEATKERNPARAMINIT, HEATKERNEXPANDPARAM, SCG, CONJGRAD


%	Copyright (c) 2010 Mauricio A. Alvarez


if kern.includeIC
    params = [kern.decay kern.diffusion kern.inverseWidthTime kern.inverseWidthSpace kern.sensitivity ...
        kern.inverseWidthSpaceIC kern.sensitivityIC];
    if nargout > 1
        names = {'decay', 'diffusion rate', 'inverse width time', 'inverse width space.', 'sensitivity', ...
            'inverse width space IC.', 'sensitivity IC'};
    end
else
    params = [kern.decay kern.diffusion kern.inverseWidthTime kern.inverseWidthSpace kern.sensitivity];
    if nargout > 1
        names = {'decay', 'diffusion rate', 'inverse width time', 'inverse width space.', 'sensitivity'};
    end
end

% if kern.includeIndSens
%     params = [params kern.sensitivitySpace];
%     if nargout > 1
%         names{5} = 'sensitivity time';
%         names{end+1} = 'sensitivity space';        
%     end
% end