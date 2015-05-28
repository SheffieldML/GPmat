function [params,names] = lmcKernExtractParam(kern)

% LMCKERNEXTRACTPARAM Extract parameters from the LMC kernel struc.
%
%	Description:
%
%	PARAM = LMCKERNEXTRACTPARAM(KERN) Extract parameters from the linear
%	model of coregionalization kernel structure into a vector of
%	parameters for optimisation.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. The vector
%	   of transforms is assumed to be empty here, any transormation of
%	   parameters is assumed to be done in the component kernels.
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%
%	[PARAM, NAMES] = LMCKERNEXTRACTPARAM(KERN) the same that before but
%	also returns the names of the parameters.
%	 Returns:
%	  PARAM - vector of parameters extracted from the kernel. The vector
%	   of transforms is assumed to be empty here, any transormation of
%	   parameters is assumed to be done in the component kernels.
%	  NAMES - cell array of strings containing parameter names.
%	 Arguments:
%	  KERN - the kernel structure containing the parameters to be
%	   extracted.
%	
%
%	See also
%	% SEEALSO MULTIKERNEXTRACTPARAM


%	Copyright (c) 2010 Mauricio A. Alvarez


% First extract the parameters of the basic kernel

fhandle = str2func([kern.basicKernelType 'KernExtractParam']);

if nargout > 1
  [paramsTemp, namesTemp] = fhandle(kern);
  namesB = cell(kern.nout, kern.rankCorregMatrix);
  for i = 1:kern.nout
      for j =1:kern.rankCorregMatrix
          namesB{i,j} = ['A(' num2str(i) ',' num2str(j) ')'];
      end
  end
  names = {namesTemp{1:kern.nParamsBK}, namesB{:}};
else
  paramsTemp = fhandle(kern);
end

% Add the parameters of the corregionalization matrix 

params = [paramsTemp(1:kern.nParamsBK)  kern.A(:)'];