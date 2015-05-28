function [paramtransformsettings, names] = ndsimKernExtractParamTransformSettings(kern)

% NDSIMKERNEXTRACTPARAMTRANSFORMSETTINGS Extract parameter transform settings from the SIM kernel structure.
%
%	Description:
%
%	PARAM = NDSIMKERNEXTRACTPARAMTRANSFORMSETTINGS(KERN) Extract
%	parameters from the single input motif kernel structure into a
%	vector of parameters for optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. If the
%	   field 'transforms' is not empty in the kernel matrix, the
%	   parameters will be transformed before optimisation (for example
%	   positive only parameters could be logged before being returned).
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = NDSIMKERNEXTRACTPARAMTRANSFORMSETTINGS(KERN)
%	Extract parameters and their names from the single input motif
%	kernel structure.
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
%
%	See also
%	% SEEALSO SIMKERNPARAMINIT, SIMKERNEXPANDPARAM, KERNEXTRACTPARAM, SCG, CONJGRAD


%	Copyright (c) 2006, 2009 Neil D. Lawrence
%	Copyright (c) 2011 Jaakko Peltonen


if isfield(kern, 'options') ...
      && isfield(kern.options, 'isNegativeS') ...
      && kern.options.isNegativeS,  
  if isfield(kern, 'gaussianInitial') && kern.gaussianInitial,
    % positive transformations for: inversewidth, initial variance
    paramtransformsettings = {kern.transforms(1).transformsettings, kern.transforms(2).transformsettings};
    if nargout > 1, names = {'inverse width', 'initial variance'}; end;
  else
    % positive transformations for: inversewidth
    paramtransformsettings = {kern.transforms(1).transformsettings};
    if nargout > 1, names = {'inverse width'}; end;
  end;
    
else
  if isfield(kern, 'gaussianInitial') && kern.gaussianInitial,
    % positive transformations for: inversewidth, kernel-variance,
    % initial variance
    paramtransformsettings = {kern.transforms(1).transformsettings, kern.transforms(2).transformsettings, kern.transforms(3).transformsettings};
    if nargout > 1, names = {'inverse width', 'variance', 'initial variance'}; end;
  else
    % positive transformations for: inversewidth, kernel-variance
    paramtransformsettings = {kern.transforms(1).transformsettings, kern.transforms(2).transformsettings};
    if nargout > 1, names = {'inverse width', 'variance'}; end;
  end;
end;


%fprintf(1, 'simKern parameters physical values:\n');
%params
%names
