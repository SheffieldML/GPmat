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
%
% MODIFICATIONS : Mauricio A. Alvarez, 2009

% KERN

if kern.isArd
    if kern.nIndFunct == 1,
        kern.precisionT =  params(1:end-1)';
    else
        kern.precisionT = reshape(params(1:end-1), size(kern.precisionT));
    end
else
    kern.precisionT =  params(1:end-1);
end
kern.sigma2Noise = params(end);
