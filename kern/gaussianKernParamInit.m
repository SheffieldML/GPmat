function kern = gaussianKernParamInit(kern, isArd)

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
%
% MODIFICATIONS : Mauricio A. Alvarez, 2009, 2010

% KERN
  
if kern.inputDimension == 0
   kern.inputDimension = 1; 
end

if isfield(kern, 'options') && isfield(kern.options, 'isArd')
    kern.isArd = kern.options.isArd;
else
    switch nargin
        case 1
            kern.isArd = false;
        case 2
            kern.isArd = isArd;
    end
end
kern.sigma2Latent = 1;
if kern.isArd
    kern.precisionU   = ones(kern.inputDimension,1);
    kern.nParams = kern.inputDimension + 1;
else
    kern.precisionU   = 1;
    kern.nParams = 2;
end
% Constrains parameters positive for optimisation.
% The variances of P need to be positive and we constrain the sensitivity
% to be positive as well
kern.transforms.index =1:kern.nParams;
kern.transforms.type = optimiDefaultConstraint('positive');
kern.isStationary = true;
