function g = multiKernDiagGradientBlock(kern, X, covGrad, i)

% MULTIKERNDIAGGRADIENTBLOCK
%
%	Description:
%
%	G = MULTIKERNDIAGGRADIENTBLOCK(KERN, X, COVGRAD, I) computes the
%	gradient with respect to parameters for the diagonal of a block of a
%	multi-output kernel given a matrix of input.
%	 Returns:
%	  G - the gradient of the kernel parameters from the kernel in the
%	   order provided by the relevant kernExtractParam commands.
%	 Arguments:
%	  KERN - the structure containing the kernel.
%	  X - first set of kernel inputs.
%	  COVGRAD - Gradient of the objective function with respect to the
%	   relevant portion of the kernel matrix.
%	  I - index of the block of the kernel to be computed.
%	
%
%	See also
%	MULTIKERNCREATE, MULTIKERNGRADIENT, MULTIKERNCOMPUTEBLOCK


%	Copyright (c) 2010 Mauricio A. Álvarez



fhandle = [kern.comp{i}.type 'KernDiagGradient'];
arg{1} = kern.comp{i};
factors = kernFactors(kern.comp{i}, 'gradfact');
fhandle = str2func(fhandle);
arg{end+1} = X;
g = fhandle(arg{:}, covGrad);
g(factors.index) = g(factors.index).*factors.val;

