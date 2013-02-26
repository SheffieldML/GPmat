function g = rbfKernDiagGradient(kern, x, covDiag)

% RBFKERNDIAGGRADIENT Compute the gradient of the RBF kernel's diagonal wrt parameters.
%
%	Description:
%
%	G = RBFKERNDIAGGRADIENT(KERN, X, FACTORS) computes the gradient of
%	functions of the diagonal of the radial basis function kernel matrix
%	with respect to the parameters of the kernel. The parameters'
%	gradients are returned in the order given by the rbfKernExtractParam
%	command.
%	 Returns:
%	  G - gradients of the relevant function with respect to each of the
%	   parameters. Ordering should match the ordering given in
%	   rbfKernExtractParam.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are computed.
%	  X - the input data for which the gradient is being computed.
%	  FACTORS - partial derivatives of the function of interest with
%	   respect to the diagonal elements of the kernel.
%	
%	
%
%	See also
%	RBFKERNPARAMINIT, KERNDIAGGRADIENT, RBFKERNEXTRACTPARAM, RBFKERNGRADIENT


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence


%	With modifications by David Luengo 2009



g = zeros(1, kern.nParams);
if isfield(kern, 'isNormalised') && (kern.isNormalised == true)
    g(1) = sqrt(kern.inverseWidth/(2*pi)) * sum(covDiag);
    g(2) = 0.5 * kern.variance / sqrt(2*pi*kern.inverseWidth) * sum(covDiag);
else
    g(1) = 0;
    g(2) = sum(covDiag);
end
