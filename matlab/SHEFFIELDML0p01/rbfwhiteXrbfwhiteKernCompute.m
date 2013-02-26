function K = rbfwhiteXrbfwhiteKernCompute(rbfKern1, rbfKern2, t1, t2)

% RBFWHITEXRBFWHITEKERNCOMPUTE Compute a cross kernel between two RBF-WHITE
%
%	Description:
%	kernels.
%
%	K = RBFWHITEXRBFWHITEKERNCOMPUTE(LFMKERN1, LFMKERN2, T1) computes
%	cross kernel terms between two RBF-WHITE kernels for the multiple
%	output kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  LFMKERN1 - the kernel structure associated with the first
%	   RBF-WHITE kernel.
%	  LFMKERN2 - the kernel structure associated with the second
%	   RBF-WHITE kernel.
%	  T1 - inputs for which kernel is to be computed.
%
%	K = RBFWHITEXRBFWHITEKERNCOMPUTE(RBFKERN1, RBFKERN2, T1, T2)
%	computes cross kernel terms between two RBF-WHITE kernels for the
%	multiple output kernel.
%	 Returns:
%	  K - block of values from kernel matrix.
%	 Arguments:
%	  RBFKERN1 - the kernel structure associated with the first
%	   RBF-WHITE kernel.
%	  RBFKERN2 - the kernel structure associated with the second
%	   RBF-WHITE kernel.
%	  T1 - row inputs for which kernel is to be computed.
%	  T2 - column inputs for which kernel is to be computed.
%	
%
%	See also
%	MULTIKERNPARAMINIT, MULTIKERNCOMPUTE, RBFWHITEKERNPARAMINIT


%	Copyright (c) 2009 David Luengo


if nargin < 4
  t2 = t1;
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end
if rbfKern1.variance ~= rbfKern2.variance
  error('Kernels cannot be cross combined if they have different variances.')
end
if rbfKern1.isStationary ~= rbfKern2.isStationary
  error('Stationary and non-stationary kernels cannot be cross combined.')
end

isStationary = rbfKern1.isStationary;

T1 = repmat(t1, 1, size(t2, 1));
T2 = repmat(t2.', size(t1, 1), 1);
deltaT = T1-T2;
indT = double(deltaT >= 0);

variance = rbfKern1.variance;
inverseWidth1 = rbfKern1.inverseWidth;
inverseWidth2 = rbfKern2.inverseWidth;
statInvWidth = inverseWidth2 + (inverseWidth1-inverseWidth2)*indT;

c = variance / sqrt(8*pi);
K = 1 - erf(sqrt(0.5*inverseWidth1*inverseWidth2) * statInvWidth ...
    .* abs(deltaT) / (inverseWidth1+inverseWidth2));
if (isStationary == false)
    K = K + erf(sqrt(0.5*inverseWidth1*inverseWidth2) * (inverseWidth1*T1 ...
        + inverseWidth2*T2) / (inverseWidth1+inverseWidth2)) - 1;
end
K = c * K ...
    .* exp((-0.5*inverseWidth1*inverseWidth2/(inverseWidth1+inverseWidth2))*(deltaT.^2));
