function g = whiteblockKernDiagGradient(kern, x, covDiag)

% WHITEBLOCKKERNDIAGGRADIENT WHITEBLOCK kernel's diagonal gradient wrt par.
%
%	Description:
%
%	G = WHITEBLOCKKERNDIAGGRADIENT(KERN, X, COVDIAG) computes the
%	gradient of the diagonal of the white noise kernel block matrix with
%	respect to the parameters of the kernel. The parameters' gradients
%	are returned in the order given by the whiteblockKernExtractParam
%	command.
%	 Returns:
%	  G - gradients of the relevant function with respect to each of the
%	   parameters. Ordering should match the ordering given in
%	   whiteblockKernExtractParam.
%	 Arguments:
%	  KERN - the kernel structure for which the gradients are computed.
%	  X - the input data for which the gradient is being computed.
%	  COVDIAG - partial derivatives of the function of interest with
%	   respect to the diagonal elements of the kernel.
%	
%
%	See also
%	WHITEBLOCKKERNPARAMINIT, KERNDIAGGRADIENT


%	Copyright (c) 2010 Mauricio A. Alvarez


g = zeros(1, kern.nout);
startOne = 1;
endOne = 0;
for i=1:kern.nout
    endOne = endOne + size(x,1);
    g(1, i) = sum(covDiag(startOne:endOne));
    startOne = endOne + 1;
end
