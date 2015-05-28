function gX = gibbsperiodicKernGradX(kern, X, X2)

% GIBBSPERIODICKERNGRADX Gradient of GIBBSPERIODIC kernel with respect to a point x.
% FORMAT
% DESC computes the gradient of the Gibbs-kernel derived periodic
% kernel with respect to the input positions. 
% ARG kern : kernel structure for which gradients are being
% computed.
% ARG x : locations against which gradients are being computed.
% RETURN g : the returned gradients. The gradients are returned in
% a matrix which is numData x numInputs x numData. Where numData is
% the number of data points and numInputs is the number of input
% dimensions in X.
%
% FORMAT
% DESC computes the gradident of the Gibbs-kernel derived periodic
% kernel with respect to the input positions where both the row
% positions and column positions are provided separately.
% ARG kern : kernel structure for which gradients are being
% computed.
% ARG x1 : row locations against which gradients are being computed.
% ARG x2 : column locations against which gradients are being computed.
% RETURN g : the returned gradients. The gradients are returned in
% a matrix which is numData2 x numInputs x numData1. Where numData1 is
% the number of data points in X1, numData2 is the number of data
% points in X2 and numInputs is the number of input
% dimensions in X.
%
% SEEALSO gibbsperiodicKernParamInit, kernGradX, gibbsperiodicKernDiagGradX
%
% COPYRIGHT : Neil D. Lawrence, 2007

% KERN


fhandle = str2func([kern.lengthScaleTransform, 'Transform']);
l = fhandle(modelOut(kern.lengthScaleFunc, X), 'atox');
l2 = fhandle(modelOut(kern.lengthScaleFunc, X2), 'atox');
gl = modelOutputGradX(kern.lengthScaleFunc, X);
gl = gl.*fhandle(repmat(l, 1, size(gl, 2)), 'gradfact');

gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
for i = 1:size(X, 1);
  gX(:, :, i) = gibbsKernGradXpoint(kern, X(i, :), X2, l(i), l2, ...
                                    gl(i, :));
end
  

function gX = gibbsKernGradXpoint(kern, x, X2, l, l2, gl)

% GIBBSKERNGRADXPOINT Gradient with respect to one point of x.


gX = zeros(size(X2));
arg = 0.5*(X2 - x);
sinarg = sin(arg);
n2 = 4*sinarg.*sinarg;
w2 = (l*l + l2.*l2);
k = kern.variance*(2*(l*l2)./w2).*exp(-n2./w2);
gX = 4*cos(arg).*sinarg./w2.*k;

base2 = k.*((l2.^4 - l^4) + 2*l*l*n2)./(w2.*w2.*l);
gX = gX + base2*gl;
