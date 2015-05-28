function kern = gaussianKernExpandParam(kern, params)

% GAUSSIANKERNEXPANDPARAM Create kernel structure from gaussian kernel's parameters.
% FORMAT
% DESC returns a gaussian kernel structure filled with the parameters in the given
%	vector. This is used as a helper function to enable parameters to be
%	optimised in, for example, the NETLAB optimisation functions.
% RETURN kern : kernel structure with the given parameters in the relevant
%	   locations.
% ARG kern : the kernel structure in which the parameters are to be
% ARG param : vector of parameters which are to be placed in the kernel
%	   structure.
%	
% SEEALSO : gaussianKernParamInit, gaussianKernExtractParam, kernExpandParam
%
% COPYRIGHT : Mauricio Alvarez and Neil D. Lawrence, 2008
%
% MODIFICATIONS : Mauricio A. Alvarez, 2009

% KERN
  
kern.sigma2Latent = params(end);
kern.precisionU =  params(1:end-1)';
