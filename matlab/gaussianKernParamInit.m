function kern = gaussianKernParamInit(kern)

% GAUSSIANKERNPARAMINIT Gaussian kernel parameter initialisation.
% The gaussian kernel used here follows the shape of a gaussian
%	distribution 
%	
%	k(x_i, x_j) =  sigma2*exp(- 0.5*(x_i - x_j)'P(x_i - x_j))
%	
%	In the above equation, P is the precision matrix and sigma2 is a variance factor. P is a diagonal matrix. 	
%
%  FORMAT
% DESC  initialises the gaussian
%	kernel structure with some default parameters.
% RETURN kern : the kernel structure with the default parameters placed in.
% ARG kern : the kernel structure which requires initialisation.
%	
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008

% KERN
  
if kern.inputDimension == 0
   kern.inputDimension = 1; 
end
kern.sigma2_u = 1;
kern.precision_u = 100*ones(kern.inputDimension,1);
kern.nParams = kern.inputDimension + 1;
% Constrains parameters positive for optimisation.
% The variances of P need to be positive and we constrain the sensitivity
% to be positive as well
kern.transforms.index =1:kern.nParams;
kern.transforms.type = optimiDefaultConstraint('positive');
kern.isStationary = true;
kern.isNormalised = false;