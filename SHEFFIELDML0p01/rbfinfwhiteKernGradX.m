function gT = rbfinfwhiteKernGradX(kern, t1, t2)

% RBFINFWHITEKERNGRADX Gradient of RBF-WHITE kernel (with integration limits
%
%	Description:
%	between minus infinity and infinity) with respect to a point t.
%
%	GT = RBFINFWHITEKERNGRADX(KERN, T1) computes the gradient of the
%	RBF-WHITE kernel with respect to the input positions.
%	 Returns:
%	  GT - the returned gradients. The gradients are returned in a
%	   matrix which is numData x numInputs x numData. Where numData is
%	   the number of data points and numInputs is the number of input
%	   dimensions in t1 (currently always one).
%	 Arguments:
%	  KERN - kernel structure for which gradients are being computed.
%	  T1 - locations against which gradients are being computed.
%
%	GT = RBFINFWHITEKERNGRADX(KERN, T1, T2) computes the gradient of the
%	RBF-WHITE kernel with respect to the input positions where both the
%	row positions and column positions are provided separately.
%	 Returns:
%	  GT - the returned gradients. The gradients are returned in a
%	   matrix which is numData2 x numInputs x numData1. Where numData1 is
%	   the number of data points in t1, numData2 is the number of data
%	   points in t2 and numInputs is the number of input dimensions in t1
%	   and t2 (currently always one).
%	 Arguments:
%	  KERN - kernel structure for which gradients are being computed.
%	  T1 - row locations against which gradients are being computed.
%	  T2 - column locations against which gradients are being computed.
%	
%
%	See also
%	% SEEALSO RBFINFWHITEKERNPARAMINIT, KERNGRADX, RBFINFWHITEKERNDIAGGRADX


%	Copyright (c) 2009 David Luengo



if nargin < 3
  t2 = t1;
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end

gT = zeros(size(t1, 1), 1, size(t2, 1));

% Parameters of the kernel required in the computation
variance = kern.variance;
invWidth = kern.inverseWidth;

for i = size(t1, 1)
    gT(i, 1, :) = - ( variance * invWidth * (t1(i)-t2) * sqrt(invWidth/pi) / 4) ...
        .* exp(- invWidth * ((t1(i)-t2).^2) / 4);
end
