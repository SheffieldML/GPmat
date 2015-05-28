function [params, names] = lfmvKernExtractParam(kern)

% LFMVKERNEXTRACTPARAM Extract parameters from the LFMV kernel structure.
%
%	Description:
%
%	PARAM = LFMVKERNEXTRACTPARAM(KERN) Extract parameters from the LFMV
%	kernel structure into a vector of parameters for optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel.
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = LFMVKERNEXTRACTPARAM(KERN) Extract parameters and
%	their names from the single input motif kernel structure.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel.
%	  NAMES - cell array of strings containing parameter names.
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%	
%
%	See also
%	% SEEALSO LFMKERNPARAMINIT, LFMKERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2010 Mauricio A. Alvarez


[params, names] = lfmKernExtractParam(kern);