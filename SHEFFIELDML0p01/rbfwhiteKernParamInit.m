function kern = rbfwhiteKernParamInit(kern)

% RBFWHITEKERNPARAMINIT RBF-WHITE kernel parameter initialisation. The RBF-
%
%	Description:
%	WHITE is a convolutional kernel defined from the convolution of a white
%	input process with a smoothing kernel given by an RBF.
%	
%	The parameters are sigma2, the process variance (kern.variance),
%	and gamma, the inverse width (kern.inverseWidth). The inverse
%	width controls how wide the basis functions are, the larger
%	gamma, the smaller the basis functions are.
%	
%
%	KERN = RBFWHITEKERNPARAMINIT(KERN) initialises the RBF-WHITE kernel
%	structure with some default parameters.
%	 Returns:
%	  KERN - the kernel structure with the default parameters placed in.
%	 Arguments:
%	  KERN - the kernel structure which requires initialisation.
%	
%	
%
%	See also
%	KERNCREATE, KERNPARAMINIT


%	Copyright (c) 2009 David Luengo
%	Copyright (c) 2009 Neil D. Lawrence


% Neil: this is odd as this kernel should be valid for any input.
if kern.inputDimension > 1
  error('RBF-WHITE kernel only valid for one-D input.')
end

kern.nParams = 2;
kern.inverseWidth = 1;
kern.variance = 1;

% Constrains parameters to be positive for optimisation.
kern.transforms.index = [1 2];
kern.transforms.type = optimiDefaultConstraint('positive');
kern.isStationary = false;
kern.positiveTime = true;
