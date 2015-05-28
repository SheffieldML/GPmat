function kern = whiteKernParamInit(kern)

% WHITEKERNPARAMINIT WHITE kernel parameter initialisation.
%
%	Description:
%	The white noise kernel arises from assuming independent Gaussian
%	noise for each point in the function. The variance of the noise is
%	given by the kern.variance parameter. The variance parameter is
%	constrained to be positive, either by an exponential
%	transformation (default) or, if the flag use_sigmoidab is set,
%	by a sigmoid transformation with a customizable output range.
%	
%	This kernel is not intended to be used independently, it is provided
%	so that it may be combined with other kernels in a compound kernel.
%	
%	
%
%	KERN = WHITEKERNPARAMINIT(KERN) initialises the white noise kernel
%	structure with some default parameters.
%	 Returns:
%	  KERN - the kernel structure with the default parameters placed in.
%	 Arguments:
%	  KERN - the kernel structure which requires initialisation.
%	
%	
%	
%
%	See also
%	CMPNDKERNPARAMINIT, KERNCREATE, KERNPARAMINIT


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence
%	Copyright (c) 2011 Jaakko Peltonen
%	Copyright (c) 2012 Antti Honkela



kern.variance = exp(-2);
kern.nParams = 1;


% The white-noise kernel can be set to use a ranged sigmoid
% (sigmoidab) transformation for the variance, instead of a plain
% exponential transformation.
if (isfield(kern,'options')) && ...
   (isfield(kern.options,'use_sigmoidab')) && ...
   (kern.options.use_sigmoidab==1),
  
  kern.transforms.index = 1;
  kern.transforms.type = 'sigmoidab'; %optimiDefaultConstraint('positive');
  kern.transforms.transformsettings = [0 1e6];
elseif (isfield(kern,'options')) && ...
   (isfield(kern.options,'boundedParam')) && kern.options.boundedParam,
  kern.transforms.index = 1;
  kern.transforms.type = optimiDefaultConstraint('bounded');
  kern.transforms.transformsettings = [0 1e6];
else
  kern.transforms.index = 1;
  kern.transforms.type = optimiDefaultConstraint('positive');
end;  
  

kern.isStationary = true;
