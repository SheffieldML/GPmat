function kern = gaussianwhiteKernParamInit(kern, isArd, nInd)

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
% ARG kern  : the kernel structure which requires initialisation.
% ARG isArd : specifies if the kernel is ARD
% ARG nInd  : number of inducing functions  	
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Mauricio A. Alvarez and Neil D. Lawrence, 2008
% 
% MODIFICATIONS : Mauricio A. Alvarez, 2009.

% KERN

% By default it assumes the kernel is ARD and only have one inducing kernel
switch nargin
    case 1
        kern.isArd = false;
        kern.nIndFunct = 1;
    case 2
        kern.isArd = isArd;
        kern.nIndFunct = 1;        
    case 3
        kern.isArd = isArd;
        kern.nIndFunct = nInd;        
    otherwise
        error('Number of inputs is incorrect')
end

kern.sigma2Noise = 1;

if kern.isArd
    if kern.nIndFunct == 1
        kern.precisionT = ones(kern.inputDimension,1);
        kern.nParams = kern.inputDimension + 1;
    else
        kern.precisionT = ones(kern.inputDimension,kern.nIndFunct);
        kern.nParams = numel(kern.precisionT) + 1;        
    end
else
    if kern.nIndFunct == 1
        kern.precisionT = 1;  % Number of rows should equal 1 if not ARD
        kern.nParams = 2;        
    else
        kern.precisionT = ones(1,kern.nIndFunct);  % Number of rows should equal 1 if not ARD
        kern.nParams = kern.nIndFunct + 1;
    end    
end
kern.transforms.index =1:kern.nParams;
kern.transforms.type = optimiDefaultConstraint('positive');
kern.isStationary = true;
