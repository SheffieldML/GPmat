function kern = pathKernParamInit(kern)

% PATHKERNPARAMINIT PATH kernel parameter initialisation.
% FORMAT
% DESC initialises the path
%  kernel structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Andrea Baisero, Carl Henrik Ek, 2013

% SHEFFIELDML


kern.cd=1.1/3;
kern.chv=.9/3;
kern.nParams=2+kern.gkern.nParams;

if(~isfield(kern,'wmatDim'))
    kern.wmatDim=0;
    kern.wmat=[];
    kern.dwdmat=[];
    kern.dwhvmat=[];
end

kern.isStationary = true;
kern.isNormalised = false;
kern.inputDimension = NaN;