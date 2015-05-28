function kern = ggwhiteKernParamInit(kern, isArd)

% GGWHITEKERNPARAMINIT GG WHITE kernel parameter initialisation.
% FORMAT
% The Gaussian Gaussian (gg) white kernel corresponds to the covariance
% matrix of an output process y(x), which is generated through a convolution
% between a latent process (white noise) with variance s^2_r and a Gaussian
% kernel G(x-s):
%	
%	y_n(x) =  sum_k int G_{nk}(x-s)u_k(s)ds
%	
% where K_{nk}(x-s) is a Gaussian kernel with precision matrix precisionG,
% and u_k(s) is an inducing function represented as a white noise process
% with variance s^2_r With this assumptions, y_n(x) is a Gaussian process
% with covariance provided by the Gaussian Gaussian white kernel.

% The kernel is designed to interoperate with the multiple output block
% kernel so that u_k(s) can be inferred given several different
% instantiations of y_n(x). precisionG is considered as diagonal.
% 
% DESC initialises the Gaussian Gaussian white kernel structure with some default parameters.
% RETURN kern : the kernel structure with the default parameters placed in.
% ARG kern : the kernel structure which requires initialisation.
%	
% SEEALSO : kernCreate, kernParamInit, ggwhiteKernCompute
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008
%
% MODIFICATIONS : Mauricio A Alvarez, 2009

% KERN

switch nargin
    case 1
        kern.isArd = false;        
    case 2
        kern.isArd = isArd;
end

if kern.isArd
    kern.precisionG = ones(kern.inputDimension,1);
    kern.sigma2Noise = 1;  % Also called variance latent or variance of noise
    kern.variance = 1;     % Also called variance output or sensitivity
    kern.nParams =kern.inputDimension + 2 ;    
    kern.transforms.index = 1:kern.inputDimension + 1;
else
    kern.precisionG  =  1;
    kern.sigma2Noise =  1;  % Also called variance latent or variance of noise
    kern.variance    =  1;     % Also called variance output or sensitivity
    kern.nParams     =  3;        
    kern.transforms.index = [1 2];
end
% The variances must be positive. As well as the sensitivity of the latent
% function.
kern.transforms.type = optimiDefaultConstraint('positive');
kern.isStationary = true;
