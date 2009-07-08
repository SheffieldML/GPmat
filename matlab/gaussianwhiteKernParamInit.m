function kern = gaussianwhiteKernParamInit(kern, nInd)

% GAUSSIANWHITEKERNPARAMINIT Gaussian white kernel parameter initialisation.
% The gaussian white kernel used here corresponds to the covariance of an 
% output function which has been obtained like the output of the
% convolution operation between a smoothing kernel function T, which follows a 
% a Gaussian form and white noise with variance s^2_r. It is given as
%  
%    k(x_i, x_j) = s^2_r \int T_i(x_i - z)T_i(x_j - z) dz   	
%	             = s^2_r N(x_i;x_j, L_{i,r}^{-1} + L_{j,r}^{-1})
%	
% where L_{i,r} corresponds to the inverse width in each direction of each
% pseudo point.
%
% FORMAT
% DESC  initialises the gaussian white kernel structure with some default 
%       parameters.
% RETURN kern : the kernel structure with the default parameters placed in.
% ARG kern : the kernel structure which requires initialisation.
% ARG nInd : number of inducing functions  	
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008
% 
% MODIFICATIONS : Mauricio A. Alvarez, 2009.

% KERN

% If isArd is true, it assumes the number of inducing functions is one. If
% isArd is false, it assumes the number of inducing functions is given by nIndFunct 

kern.isArd = true;
if nargin < 2
    kern.nIndFunct = 20;              % Number of inducing functions.
else
    kern.nIndFunct = nInd;           % Number of inducing functions.
end
kern.sigma2Noise = 1;
if kern.isArd
    kern.precisionT = ones(kern.inputDimension,1);
    kern.nParams = kern.inputDimension + 1;
else
    kern.precisionT = 1e-2*ones(1,kern.nIndFunct);  % Number of rows should equal 1 if not ARD
    kern.nParams = kern.nIndFunct + 1;    
end
kern.transforms.index =1:kern.nParams;
kern.transforms.type = optimiDefaultConstraint('positive');
kern.isStationary = true;
