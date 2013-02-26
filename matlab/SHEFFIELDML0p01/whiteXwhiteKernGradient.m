function [g1, g2] = whiteXwhiteKernGradient(whiteKern1, whiteKern2, x1, x2, covGrad)

% WHITEXWHITEKERNGRADIENT Compute a cross gradient between two WHITE kernels.
%
%	Description:
%
%	[G1, G2] = WHITEXWHITEKERNGRADIENT(WHITEKERN1, WHITEKERN2, X,
%	COVGRAD) computes cross gradient of parameters of a cross kernel
%	between two white noise kernels for the multiple output kernel.
%	 Returns:
%	  G1 - gradient of the parameters of the first kernel, for ordering
%	   see whiteKernExtractParam.
%	  G2 - gradient of the parameters of the second kernel, for ordering
%	   see whiteKernExtractParam.
%	 Arguments:
%	  WHITEKERN1 - the kernel structure associated with the first white
%	   noise kernel.
%	  WHITEKERN2 - the kernel structure associated with the second white
%	   noise kernel.
%	  X - inputs for which kernel is to be computed.
%	  COVGRAD - gradient of the objective function with respect to the
%	   elements of the cross kernel matrix.
%
%	[G1, G2] = WHITEXWHITEKERNGRADIENT(WHITEKERN1, WHITEKERN2, X1, X2,
%	COVGRAD) computes cross kernel terms between two white noise kernels
%	for the multiple output kernel.
%	 Returns:
%	  G1 - gradient of the parameters of the first kernel, for ordering
%	   see whiteKernExtractParam.
%	  G2 - gradient of the parameters of the second kernel, for ordering
%	   see whiteKernExtractParam.
%	 Arguments:
%	  WHITEKERN1 - the kernel structure associated with the first white
%	   noise kernel.
%	  WHITEKERN2 - the kernel structure associated with the second white
%	   noise kernel.
%	  X1 - row inputs for which kernel is to be computed.
%	  X2 - column inputs for which kernel is to be computed.
%	  COVGRAD - gradient of the objective function with respect to the
%	   elements of the cross kernel matrix.
%	
%
%	See also
%	MULTIKERNPARAMINIT, MULTIKERNCOMPUTE, WHITEKERNPARAMINIT, WHITEKERNEXTRACTPARAM


%	Copyright (c) 2006 Neil D. Lawrence


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