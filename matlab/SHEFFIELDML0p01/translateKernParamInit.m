function kern = translateKernParamInit(kern)

% TRANSLATEKERNPARAMINIT TRANSLATE kernel parameter initialisation.
%
%	Description:
%	This kernel allows you to take any other kernel in the tool box and
%	translate with respect to the input space. With stationary kernels this is
%	obviously pointless, but for non-stationary kernels (such as the LINEAR
%	kernel or the MLP kernel) it moves the non-stationarities around the input
%	space which is potentially useful.
%	
%	
%
%	KERN = TRANSLATEKERNPARAMINIT(KERN) initialises the input space
%	translation kernel structure with some default parameters.
%	 Returns:
%	  KERN - the kernel structure with the default parameters placed in.
%	 Arguments:
%	  KERN - the kernel structure which requires initialisation.
%	
%
%	See also
%	CMPNDKERNPARAMINIT, KERNCREATE, KERNPARAMINIT


%	Copyright (c) 2007 Neil D. Lawrence


kern = cmpndKernParamInit(kern);

% Add in the translation parameters
kern.nParams = kern.nParams + kern.inputDimension;
kern.centre = zeros(1, kern.inputDimension);
