function kern = gaussianwhiteKernExpandParam(kern, params)

% GAUSSIANWHITEKERNEXPANDPARAM Create kernel structure from gaussian white 
%                              kernel's parameters.
% FORMAT
% DESC returns a gaussian white kernel structure filled with the parameters in the given
%	vector. This is used as a helper function to enable parameters to be
%	optimised in, for example, the NETLAB optimisation functions.
% RETURN kern : kernel structure with the given parameters in the relevant
%	   locations.
% ARG kern : the kernel structure in which the parameters are to be
% ARG param : vector of parameters which are to be placed in the kernel
%	   structure.
%	
% SEEALSO : gaussianwhiteKernParamInit, gaussianwhiteKernExtractParam, kernExpandParam
%
% COPYRIGHT : Mauricio Alvarez and Neil D. Lawrence, 2008

% KERN
  
kern.sigma2Noise = params(end);
kern.precisionT =  params(1:end-1)';