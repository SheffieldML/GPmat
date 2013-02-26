function kern = gaussianwhiteKernExpandParam(kern, params)

% GAUSSIANWHITEKERNEXPANDPARAM Create kernel structure from gaussian white
%
%	Description:
%	kernel's parameters.
%
%	KERN = GAUSSIANWHITEKERNEXPANDPARAM(KERN, PARAM) returns a gaussian
%	white kernel structure filled with the parameters in the given
%	vector. This is used as a helper function to enable parameters to be
%	optimised in, for example, the NETLAB optimisation functions.
%	 Returns:
%	  KERN - kernel structure with the given parameters in the relevant
%	   locations.
%	 Arguments:
%	  KERN - the kernel structure in which the parameters are to be
%	  PARAM - vector of parameters which are to be placed in the kernel
%	   structure.
%	
%	
%
%	See also
%	GAUSSIANWHITEKERNPARAMINIT, GAUSSIANWHITEKERNEXTRACTPARAM, KERNEXPANDPARAM


%	Copyright (c) 2008 Mauricio Alvarez and Neil D. Lawrence


%	With modifications by Mauricio A. Alvarez 2009


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