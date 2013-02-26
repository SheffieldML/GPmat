function kern = velotransKernParamInit(kern)

% VELOTRANSKERNPARAMINIT VELOTRANS kernel parameter initialisation.
%
%	Description:
%	This kernel allows you to take any other kernel in the tool box and
%	translate with respect to the input space with respect to time. The translation allows a moving field such as clouds moving over a land surface.
%	
%	
%
%	KERN = VELOTRANSKERNPARAMINIT(KERN) initialises the velocity
%	translate kernel structure with some default parameters.
%	 Returns:
%	  KERN - the kernel structure with the default parameters placed in.
%	 Arguments:
%	  KERN - the kernel structure which requires initialisation.
%	
%
%	See also
%	CMPNDKERNCREATE, KERNCREATE, KERNPARAMINIT


%	Copyright (c) 2011 Neil D. Lawrence


kern = cmpndKernParamInit(kern);

% Add in the translation parameters
kern.nParams = kern.nParams + kern.inputDimension-1;
kern.velocity = zeros(1, kern.inputDimension-1);
