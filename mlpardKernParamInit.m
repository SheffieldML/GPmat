function kern = mlpardKernParamInit(kern)

% MLPARDKERNPARAMINIT MLPARD kernel parameter initialisation.
% The automatic relevance determination variant of the multi-layer
% perceptron kernel (MLP) is sometimes also known as the arcsin
% kernel. It is derived by considering the infinite limit of a neural
% network with erf() functions as the hidden layer activations.
%
% k(x_i, x_j) = sigma2*arcsin((w*x_i'*A*x_j + b) ...
%    /sqrt((w*x_i'*A*x_i + b + 1)*(w*x_j'*A*x_j + b + 1)))
% 
% As well as the diagonal matrix of input scales (kern.inputScales:
% constrained to be between zero and one) there are three parameters,
% the kernel variance (sigma2) stored in kern.variance, the variance
% of the weights in the neural network (w) stored in
% kern.weightVariance and the variances of the biases in the neural
% network (kern.biasVariance).
% 
% SEEALSO : mlpKernParamInit
%
% FORMAT
% DESC initialises the automatic relevance determination multi-layer perceptron
%  kernel structure with some default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% KERN


kern.weightVariance = 10;
kern.biasVariance = 10;
kern.variance = 1;
% These parameters are restricted to lie between 0 and 1.
kern.inputScales = 0.999*ones(1, kern.inputDimension);
kern.nParams = 3 + kern.inputDimension;

kern.transforms(1).index = [1 2 3];
kern.transforms(1).type = optimiDefaultConstraint('positive');
kern.transforms(2).index = [4:kern.nParams];
kern.transforms(2).type = optimiDefaultConstraint('zeroone');

kern.isStationary = false;
