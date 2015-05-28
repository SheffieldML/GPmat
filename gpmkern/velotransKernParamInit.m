function kern = velotransKernParamInit(kern)

% VELOTRANSKERNPARAMINIT VELOTRANS kernel parameter initialisation.
% This kernel allows you to take any other kernel in the tool box and
% translate with respect to the input space with respect to time. The translation allows a moving field such as clouds moving over a land surface.
%
% SEEALSO : cmpndKernCreate
%
% FORMAT
% DESC initialises the velocity translate
%  kernel structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2011

% KERN

kern = cmpndKernParamInit(kern);

% Add in the translation parameters
kern.nParams = kern.nParams + kern.inputDimension-1;
kern.velocity = zeros(1, kern.inputDimension-1);
