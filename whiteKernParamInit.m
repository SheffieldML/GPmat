function kern = whiteKernParamInit(kern)

% WHITEKERNPARAMINIT WHITE kernel parameter initialisation.
% The white noise kernel arises from assuming independent Gaussian
% noise for each point in the function. The variance of the noise is
% given by the kern.variance parameter. The variance parameter is
% constrained to be positive, either by an exponential
% transformation (default) or, if the flag use_sigmoidab is set, 
% by a sigmoid transformation with a customizable output range.
% 
% This kernel is not intended to be used independently, it is provided
% so that it may be combined with other kernels in a compound kernel.
%
% SEEALSO : cmpndKernParamInit
%
% FORMAT
% DESC initialises the white noise
%  kernel structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
%
% COPYRIGHT : Jaakko Peltonen, 2011
%
% COPYRIGHT : Antti Honkela, 2012

% KERN


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
elseif (isfield(kern,'options')) && ...
      (isfield(kern.options,'paramTransform')),
  kern.transforms.index = 1;
  switch kern.options.paramTransform,
   case 'sigmoidab',
    kern.transforms.type = 'sigmoidab';
    kern.transforms.transformsettings = [0 1e6];
   case 'bounded',
    kern.transforms.type = optimiDefaultConstraint('bounded');
    kern.transforms.transformsettings = [0 1e6];
   case 'identity',
    kern.transforms.type = 'identity';
    kern.transforms.transformsettings = [0 1e6];
   case 'positive',
    kern.transforms.type = optimiDefaultConstraint('positive');
   case 'none',
   otherwise,
    error('Unknown paramTransform');
  end
else
  kern.transforms.index = 1;
  kern.transforms.type = optimiDefaultConstraint('positive');
end;  

kern.isStationary = true;
