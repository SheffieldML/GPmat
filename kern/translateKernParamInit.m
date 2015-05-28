function kern = translateKernParamInit(kern)

% TRANSLATEKERNPARAMINIT TRANSLATE kernel parameter initialisation.
% This kernel allows you to take any other kernel in the tool box and
% translate with respect to the input space. With stationary kernels this is
% obviously pointless, but for non-stationary kernels (such as the LINEAR
% kernel or the MLP kernel) it moves the non-stationarities around the input
% space which is potentially useful.
%
% SEEALSO : cmpndKernParamInit
%
% FORMAT
% DESC initialises the input space translation
%  kernel structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2007

% KERN

kern = cmpndKernParamInit(kern);

% Add in the translation parameters
kern.nParams = kern.nParams + kern.inputDimension;
kern.centre = zeros(1, kern.inputDimension);
