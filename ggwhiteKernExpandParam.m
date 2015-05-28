function kern = ggwhiteKernExpandParam(kern, params)

% GGWHITEKERNEXPANDPARAM Create kernel structure from GG white kernel's parameters.
% FORMAT
% DESC returns a gaussian gaussian white
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
% SEEALSO : ggwhiteKernParamINit, ggwhiteKernExtractParam, kernExpandParam
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008
%
% MODIFICATIONS : Mauricio A Alvarez, 2009

% KERN

kern.precisionG = params(1:end-2)';
kern.sigma2Noise = params(end-1);
kern.variance = params(end);

