function [g1, g2] = rbfwhiteXrbfwhiteKernGradient(rbfKern1, rbfKern2, t1, varargin)

% RBFWHITEXRBFWHITEKERNGRADIENT Compute a cross gradient between two
% RBF-WHITE kernels.
% FORMAT
% DESC computes cross gradient of parameters of a cross kernel
% between two RBF-WHITE kernels for the multiple output kernel. 
% ARG rbfKern1 : the kernel structure associated with the first RBF-WHITE
% kernel.
% ARG rbfKern2 : the kernel structure associated with the second RBF-WHITE
% kernel.
% ARG t1 : inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see rbfwhiteKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see rbfwhiteKernExtractParam.
%
% FORMAT
% DESC computes cross kernel terms between two RBF-WHITE kernels for
% the multiple output kernel. 
% ARG rbfKern1 : the kernel structure associated with the first RBF-WHITE
% kernel.
% ARG rbfKern2 : the kernel structure associated with the second RBF-WHITE
% kernel.
% ARG t1 : row inputs for which kernel is to be computed.
% ARG t2 : column inputs for which kernel is to be computed.
% ARG covGrad : gradient of the objective function with respect to
% the elements of the cross kernel matrix.
% RETURN g1 : gradient of the parameters of the first kernel, for
% ordering see rbfwhiteKernExtractParam.
% RETURN g2 : gradient of the parameters of the second kernel, for
% ordering see rbfwhiteKernExtractParam.
%
% SEEALSO : multiKernParamInit, multiKernCompute, rbfwhiteKernParamInit,
% rbfwhiteKernExtractParam
%
% COPYRIGHT : David Luengo, 2009

% KERN


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
