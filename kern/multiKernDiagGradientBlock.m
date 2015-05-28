function g = multiKernDiagGradientBlock(kern, X, covGrad, i)

% MULTIKERNDIAGGRADIENTBLOCK
% FORMAT
% DESC computes the gradient with respect to parameters for the diagonal of
% a block of a multi-output kernel given a matrix of input. 
% ARG kern : the structure containing the kernel.
% ARG X : first set of kernel inputs.
% ARG covGrad : Gradient of the objective function with respect to
% the relevant portion of the kernel matrix.
% ARG i : index of the block of the kernel to be computed.
% RETURN g : the gradient of the kernel parameters from the 
% kernel in the order provided by the relevant kernExtractParam commands.
%
% SEEALSO : multiKernCreate, multiKernGradient, multiKernComputeBlock
%
% COPYRIGHT : Mauricio A. Álvarez, 2010

% KERN


fhandle = [kern.comp{i}.type 'KernDiagGradient'];
arg{1} = kern.comp{i};
factors = kernFactors(kern.comp{i}, 'gradfact');
fhandle = str2func(fhandle);
arg{end+1} = X;
g = fhandle(arg{:}, covGrad);
g(factors.index) = g(factors.index).*factors.val;

