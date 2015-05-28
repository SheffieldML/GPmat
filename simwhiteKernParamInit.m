function kern = simwhiteKernParamInit(kern)

% SIMWHITEKERNPARAMINIT SIM-WHITE kernel parameter initialisation.
% FORMAT
% DESC initialises the SIM-White (Single Input Motif - White)
%  kernel structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : David Luengo, 2009

% KERN


if kern.inputDimension > 1
  error('SIM-WHITE kernel only valid for one-D input.')
end

kern.nParams = 3;
kern.decay = 1;
kern.variance = 1;
kern.sensitivity = 1;

kern.delay = 0;
kern.initVal = 1;

kern.transforms.index = [1 2];
kern.transforms.type = optimiDefaultConstraint('positive');

kern.isStationary = false;
kern.positiveTime = true;
