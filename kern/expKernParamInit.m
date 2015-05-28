function kern = expKernParamInit(kern)

% EXPKERNPARAMINIT EXP kernel parameter initialisation.
% FORMAT
% DESC initialises the exponentiated
%  kernel structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN

kern.argument = kernParamInit(kern.argument);
kern.nParams = kern.argument.nParams + 1;
kern.variance = 1;
kern.transforms.index = [1];
kern.transforms.type = optimiDefaultConstraint('positive');
if isfield(kern.argument, 'isStationary') 
  kern.isStationary = kern.argument.isStationary;
else
  kern.isStationary = false;
end
% Deal with fact that white variance is exponentiated.
if isfield(kern.argument, 'whiteVariance')
  whiteVar = kern.argument.whiteVariance;
  kern.whiteVariance = kern.variance*(exp(2*whiteVar)-exp(whiteVar));
end
