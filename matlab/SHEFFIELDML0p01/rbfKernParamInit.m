function kern = rbfKernParamInit(kern)

% RBFKERNPARAMINIT RBF kernel parameter initialisation.
%
%	Description:
%	The radial basis function kernel (RBF) is sometimes also known as
%	the squared exponential kernel. It is a very smooth non-linear
%	kernel and is a popular choice for generic use.
%	
%	k(x_i, x_j) = sigma2 * exp(-gamma/2 *(x_i - x_j)'*(x_i - x_j))
%	
%	The parameters are sigma2, the process variance (kern.variance)
%	and gamma, the inverse width (kern.inverseWidth). The inverse
%	width controls how wide the basis functions are, the larger
%	gamma, the smaller the basis functions are.
%	
%	There is also an automatic relevance determination version of
%	this kernel provided.
%	
%	
%
%	KERN = RBFKERNPARAMINIT(KERN) initialises the radial basis function
%	kernel structure with some default parameters.
%	 Returns:
%	  KERN - the kernel structure with the default parameters placed in.
%	 Arguments:
%	  KERN - the kernel structure which requires initialisation.
%	
%	
%
%	See also
%	RBFARDKERNPARAMINIT, KERNCREATE, KERNPARAMINIT


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence
%	Copyright (c) 2012 Antti Honkela



kern.inverseWidth = 1;
kern.variance = 1;
kern.nParams = 2;

% Constrains parameters positive for optimisation.
if (isfield(kern,'options')) && ...
   (isfield(kern.options,'boundedParam')) && kern.options.boundedParam,
  for k=1:2,
    kern.transforms(k).index = k;
    kern.transforms(k).type = optimiDefaultConstraint('bounded');
    kern.transforms(k).transformsettings = [0 1e6];
  end
else
  kern.transforms.index = [1 2];
  kern.transforms.type = optimiDefaultConstraint('positive');
end
kern.isStationary = true;
kern.isNormalised = false;
