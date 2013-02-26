function [g1, g2] = rbfwhiteXrbfwhiteKernGradient(rbfKern1, rbfKern2, t1, varargin)

% RBFWHITEXRBFWHITEKERNGRADIENT Compute a cross gradient between two
%
%	Description:
%	RBF-WHITE kernels.
%
%	[G1, G2] = RBFWHITEXRBFWHITEKERNGRADIENT(RBFKERN1, RBFKERN2, T1,
%	COVGRAD) computes cross gradient of parameters of a cross kernel
%	between two RBF-WHITE kernels for the multiple output kernel.
%	 Returns:
%	  G1 - gradient of the parameters of the first kernel, for ordering
%	   see rbfwhiteKernExtractParam.
%	  G2 - gradient of the parameters of the second kernel, for ordering
%	   see rbfwhiteKernExtractParam.
%	 Arguments:
%	  RBFKERN1 - the kernel structure associated with the first
%	   RBF-WHITE kernel.
%	  RBFKERN2 - the kernel structure associated with the second
%	   RBF-WHITE kernel.
%	  T1 - inputs for which kernel is to be computed.
%	  COVGRAD - gradient of the objective function with respect to the
%	   elements of the cross kernel matrix.
%
%	[G1, G2] = RBFWHITEXRBFWHITEKERNGRADIENT(RBFKERN1, RBFKERN2, T1, T2,
%	COVGRAD) computes cross kernel terms between two RBF-WHITE kernels
%	for the multiple output kernel.
%	 Returns:
%	  G1 - gradient of the parameters of the first kernel, for ordering
%	   see rbfwhiteKernExtractParam.
%	  G2 - gradient of the parameters of the second kernel, for ordering
%	   see rbfwhiteKernExtractParam.
%	 Arguments:
%	  RBFKERN1 - the kernel structure associated with the first
%	   RBF-WHITE kernel.
%	  RBFKERN2 - the kernel structure associated with the second
%	   RBF-WHITE kernel.
%	  T1 - row inputs for which kernel is to be computed.
%	  T2 - column inputs for which kernel is to be computed.
%	  COVGRAD - gradient of the objective function with respect to the
%	   elements of the cross kernel matrix.
%	rbfwhiteKernExtractParam
%	
%
%	See also
%	MULTIKERNPARAMINIT, MULTIKERNCOMPUTE, RBFWHITEKERNPARAMINIT, 


%	Copyright (c) 2009 David Luengo



if nargin < 5
    t2 = t1;
else
    t2 = varargin{1};
end
covGrad = varargin{end};

if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end
if rbfKern1.variance ~= rbfKern2.variance
  error('Kernels cannot be cross combined if they have different variances.')
end
if rbfKern1.isStationary ~= rbfKern2.isStationary
  error('Stationary and non-stationary kernels cannot be cross combined.')
end

g1 = zeros(1, rbfKern1.nParams);
g2 = zeros(1, rbfKern2.nParams);

T1 = repmat(t1, 1, size(t2, 1));
T2 = repmat(t2.', size(t1, 1), 1);
deltaT = T1 - T2;
indT = (T1 >= T2);

% Parameters required for further computations
isStationary = rbfKern1.isStationary;
variance = rbfKern1.variance;
invWidth1 = rbfKern1.inverseWidth;
invWidth2 = rbfKern2.inverseWidth;
invWidthStat = invWidth2 + (invWidth1-invWidth2)*indT;

% Computation of a normalised kernel with unitary variance
rbfKern1.variance = 1;
rbfKern2.variance = 1;
K = rbfwhiteXrbfwhiteKernCompute(rbfKern1, rbfKern2, t1, t2);

% Gradient w.r.t. the inverse widths
gradNu1 = (invWidth2/(invWidth1+invWidth2))^2;
gradNu1Erf1 = 0.5*((3*invWidth1*(invWidth2^2)+(invWidth1^2)*invWidth2)*T1 ...
    + ((invWidth2^3)-invWidth1*(invWidth2^2))*T2) ...
    / (sqrt(invWidth1*invWidth2)*(invWidth1+invWidth2)^2);
gradNu2Erf1 = 0.5 * invWidth2 * abs(deltaT) .* (invWidth2*(invWidth2-invWidth1) ...
    + (4*invWidth1*invWidth2 + (invWidth1^2) - (invWidth2^2)) * indT) ...
    / (sqrt(invWidth1*invWidth2)*(invWidth1+invWidth2)^2);
gradErf1 = sqrt(2/pi) * (gradNu1Erf1...
        .* exp(-0.5*invWidth1*invWidth2*((invWidth1*T1+invWidth2*T2)/(invWidth1+invWidth2)).^2) ...
    - exp(-0.5*invWidth1*invWidth2*(invWidthStat.*abs(deltaT)/(invWidth1+invWidth2)).^2) ...
        .* gradNu2Erf1);
g1(1) = variance * sum(sum( (-0.5 * (deltaT.^2) .* K * gradNu1 ...
    + 1/sqrt(8*pi) * gradErf1 ...
        .* exp(-0.5*invWidth1*invWidth2*(deltaT.^2)/(invWidth1+invWidth2))) ...
    .* covGrad));
gradNu2 = (invWidth1/(invWidth1+invWidth2))^2;
gradNu1Erf2 = 0.5*((3*(invWidth1^2)*invWidth2+invWidth1*(invWidth2^2))*T2 ...
    + ((invWidth1^3)-(invWidth1^2)*invWidth2)*T1) ...
    / (sqrt(invWidth1*invWidth2)*(invWidth1+invWidth2)^2);
gradNu2Erf2 = 0.5 * invWidth1 * abs(deltaT) .* (invWidth2*(3*invWidth1+invWidth2) ...
    + ((invWidth1^2) - 4*invWidth1*invWidth2 - (invWidth2^2)) * indT) ...
    / (sqrt(invWidth1*invWidth2)*(invWidth1+invWidth2)^2);
gradErf2 = sqrt(2/pi) * (gradNu1Erf2...
        .* exp(-0.5*invWidth1*invWidth2*((invWidth1*T1+invWidth2*T2)/(invWidth1+invWidth2)).^2) ...
    - exp(-0.5*invWidth1*invWidth2*(invWidthStat.*abs(deltaT)/(invWidth1+invWidth2)).^2) ...
        .* gradNu2Erf2);
g2(1) = variance * sum(sum( (-0.5 * (deltaT.^2) .* K * gradNu2 ...
    + 1/sqrt(8*pi) * gradErf2 ...
        .* exp(-0.5*invWidth1*invWidth2*(deltaT.^2)/(invWidth1+invWidth2))) ...
    .* covGrad));

% Gradient w.r.t. sigma_r^2
g1(2) = sum(sum(K .* covGrad));
g2(2) = 0; % Otherwise it is counted twice
