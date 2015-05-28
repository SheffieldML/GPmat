function kern = rbfhKernParamInit(kern)

% RBFHKERNPARAMINIT RBFH kernel parameter initialisation.
% The radial basis function for the heat kernel (RBFH) is a strictly two
% input dimension covariance function, in which each of both inputs is
% described by and RBF covariance. It is supposed to be use with the HEAT
% kernel in the latent force model approach.
% FORMAT
% DESC initialises the radial basis function for the heat kernel structure 
% with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : heatKernParamInit.m, kernCreate, kernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez, 2010

% KERN

if kern.inputDimension ~= 2
  error('RBFH kernel only valid for 2-D inputs.')
end

kern.inverseWidthTime = 1;
kern.inverseWidthSpace = 1000;

kern.rbf.inputDimension = 1;
if isfield(kern, 'options')
    kern.rbf.options = kern.options;
end
kern.rbf = rbfKernParamInit(kern.rbf);
kern.nParams = 2;
kern.transforms.type = optimiDefaultConstraint('positive');
kern.transforms.index = [1 2]; 

if isfield(kern, 'options') ...
        && isfield(kern.options, 'isStationary') ...
        && kern.options.isStationary,
   kern.isStationary = true;
else
   kern.isStationary = false;
end
