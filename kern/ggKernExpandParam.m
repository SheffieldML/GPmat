function kern = ggKernExpandParam(kern, params)

% GGKERNEXPANDPARAM Create kernel structure from GG kernel's parameters.
% FORMAT
% DESC returns a gaussian gaussian
%	kernel structure filled with the parameters in the given vector.
%	This is used as a helper function to enable parameters to be
%	optimised in, for example, the NETLAB optimisation functions.
% RETURN kern : kernel structure with the given parameters in the relevant
%	   locations.
% ARG kern : the kernel structure in which the parameters are to be
%	   placed.
% ARG param : vector of parameters which are to be placed in the kernel
%	   structure.
%	
% SEEALSO : ggKernParamINit, ggKernExtractParam, kernExpandParam
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008
%
% MODIFICATIONS : Mauricio A. Alvarez, 2009

% KERN

sizeP = size(kern.precisionU,1);
kern.precisionU = params(1:sizeP)';
kern.precisionG = params(sizeP+1:2*sizeP)';
kern.sigma2Latent = params(end-1);
kern.sensitivity = params(end);
