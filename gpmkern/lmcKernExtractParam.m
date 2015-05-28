function [params,names] = lmcKernExtractParam(kern)

% LMCKERNEXTRACTPARAM Extract parameters from the LMC kernel struc.
% FORMAT
% DESC Extract parameters from the linear model of coregionalization kernel 
% structure into a vector of parameters for optimisation. 
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the
% kernel. The vector of transforms is assumed to be empty here, any
% transormation of parameters is assumed to be done in the
% component kernels.
%
% FORMAT
% DESC the same that before but also returns the names of the parameters.
% ARG kern : the kernel structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the
% kernel. The vector of transforms is assumed to be empty here, any
% transormation of parameters is assumed to be done in the
% component kernels.
% RETURN names : cell array of strings containing parameter names.
%
% SEEALSO multiKernExtractParam
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

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
