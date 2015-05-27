function kern = biasKernParamInit(kern)

% BIASKERNPARAMINIT BIAS kernel parameter initialisation.
% The bias kernel arises from placing a prior over the bias with a
% variance given by the kern.variance parameter. It allows the
% output function to move up and down. 
%
% This kernel is not intended to be used independently, it is
% provided so that it may be combined with other kernels in a
% compound kernel.
%
% SEEALSO : cmpndKernParamInit
%
% FORMAT
% DESC initialises the bias
%  kernel structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006
%
% COPYRIGHT : Antti Honkela, 2012

% KERN


kern.variance = exp(-2);
kern.nParams = 1;

kern.transforms.index = 1;
if (isfield(kern,'options')) && ...
   (isfield(kern.options,'boundedParam')) && kern.options.boundedParam,
  kern.transforms.type = optimiDefaultConstraint('bounded');
  kern.transforms.transformsettings = [0 1e6];
else
  kern.transforms.type = optimiDefaultConstraint('positive');
end

kern.isStationary = true;
