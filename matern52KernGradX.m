function gX = matern52KernGradX(kern, X, X2)

% MATERN52KERNGRADX Gradient of MATERN52 kernel with respect to input locations.
% FORMAT
% DESC computes the gradident of the matern kernel with nu=5/2
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
% SEEALSO matern52KernParamInit, kernGradX, matern52KernDiagGradX
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN


gX = zeros(size(X2, 1), size(X2, 2), size(X, 1));
warnState = warning('query', 'MATLAB:divideByZero');
warning('off', 'MATLAB:divideByZero');
for i = 1:size(X, 1);
  gX(:, :, i) = matern52KernGradXpoint(kern, X(i, :), X2);
end
warning(warnState.state, 'MATLAB:divideByZero');
  

function gX = matern52KernGradXpoint(kern, x, X2)

% MATERN52KERNGRADXPOINT Gradient with respect to one point of x.

gX = zeros(size(X2));
n2 = dist2(X2, x);
wi2 = (5/(kern.lengthScale*kern.lengthScale));
n2wi2 = n2*wi2;
sqrtn2wi2 = sqrt(n2*wi2);
K = kern.variance*(1+sqrtn2wi2+n2wi2/3).*exp(-sqrtn2wi2);
for i = 1:size(x, 2)
  ratio = (X2(:, i) - x(i))./sqrtn2wi2;
  ratio(find(isnan(ratio)))=1;
  gX(:, i) = wi2*ratio.*(K-kern.variance.*exp(-sqrtn2wi2))-kern.variance*2*wi2/3*(X2(:, ...
                                                    i)-x(i)).*exp(-sqrtn2wi2);
end
