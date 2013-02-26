
% RBFPERIODICKERNDIAGGRADIENT Compute the gradient of the RBFPERIODIC kernel's diagonal wrt parameters.
%
%	Description:
%
%	G = RBFPERIODICKERNDIAGGRADIENT(KERN, X, FACTORS) computes the
%	gradient of functions of the diagonal of the RBF derived periodic
%	kernel matrix with respect to the parameters of the kernel. The
%	parameters' gradients are returned in the order given by the
%	rbfperiodicKernExtractParam command.
%	 Returns:
%	  G - gradients of the relevant function with respect to each of the
%	   parameters. Ordering should match the ordering given in
%	   rbfperiodicKernExtractParam.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are computed.
%	  X - the input data for which the gradient is being computed.
%	  FACTORS - partial derivatives of the function of interest with
%	   respect to the diagonal elements of the kernel.
%	
%
%	See also
%	RBFPERIODICKERNPARAMINIT, KERNDIAGGRADIENT, RBFPERIODICKERNEXTRACTPARAM, RBFPERIODICKERNGRADIENT


%	Copyright (c) 2007 Neil D. Lawrence


g = zeros(1, kern.nParams);
g(1) = 0;
g(2) = sum(covDiag);
