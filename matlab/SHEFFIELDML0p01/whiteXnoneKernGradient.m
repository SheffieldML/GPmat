function [g1, g2] = whiteXnoneKernGradient(whiteKern, noneKern, x1, x2, covGrad)

% WHITEXNONEKERNGRADIENT Compute a cross gradient between WHITE and DUMMY kernels.
%
%	Description:
%
%	[G1, G2] = WHITEXNONEKERNGRADIENT(WHITEKERN, NONEKERN, X, COVGRAD)
%	computes cross gradient of parameters of a cross kernel between a
%	white noise kernel and a dummy kernel.
%	 Returns:
%	  G1 - gradient of the parameters of the first kernel, for ordering
%	   see whiteKernExtractParam.
%	  G2 - gradient of the parameters of the second kernel, for ordering
%	   see whiteKernExtractParam.
%	 Arguments:
%	  WHITEKERN - the kernel structure associated with the white noise
%	   kernel.
%	  NONEKERN - the dummy kernel structure.
%	  X - inputs for which kernel is to be computed.
%	  COVGRAD - gradient of the objective function with respect to the
%	   elements of the cross kernel matrix.
%
%	[G1, G2] = WHITEXNONEKERNGRADIENT(WHITEKERN, NONEKERN, X1, X2,
%	COVGRAD) computes cross kernel terms between a white noise kernel
%	and a dummy kernel for the multiple output kernel.
%	 Returns:
%	  G1 - gradient of the parameters of the first kernel, for ordering
%	   see whiteKernExtractParam.
%	  G2 - gradient of the parameters of the second kernel, for ordering
%	   see whiteKernExtractParam.
%	 Arguments:
%	  WHITEKERN - the kernel structure associated with the white noise
%	   kernel.
%	  NONEKERN - the dummy kernel structure.
%	  X1 - row inputs for which kernel is to be computed.
%	  X2 - column inputs for which kernel is to be computed.
%	  COVGRAD - gradient of the objective function with respect to the
%	   elements of the cross kernel matrix.
%	
%
%	See also
%	MULTIKERNPARAMINIT, MULTIKERNCOMPUTE, WHITEKERNPARAMINIT, WHITEKERNEXTRACTPARAM


%	Copyright (c) 2008 Neil D. Lawrence


arg{1}=x1;
if nargin < 5
  covGrad = x2;
  x2 = x1;
else
  arg{2}=x2;
end
% if size(x1, 2) > 1 | size(x2, 2) > 1
%   error('Input can only have one column');
% end
g1 = 0;
g2 = 0;