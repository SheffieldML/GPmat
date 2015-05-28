function gX = gibbsKernGradX(kern, X, X2)

% GIBBSKERNGRADX Gradient of GIBBS kernel with respect to input locations.
% FORMAT
% DESC computes the gradident of the Mark Gibbs's non-stationary
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
% SEEALSO gibbsKernParamInit, kernGradX, gibbsKernDiagGradX
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN

%  pretty(simple(diff((2*(l_i*l_j)/(l_i*l_i + l_j*l_j))^(d/2)*exp(-r*r/(l_i*l_i + l_j*l_j)),l_i)))

%                                         2
%       /    l_i l_j  \(1/2 d)           r               4        4        2  2
%   1/2 |2 -----------|        exp(- -----------) (-d l_i  + d l_j  + 4 l_i  r )
%       |     2      2|                 2      2
%       \  l_i  + l_j /              l_i  + l_j

%            /      2      2 2
%           /  ((l_i  + l_j )  l_i)
%          /


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
n2 = dist2(X2, x);
w2 = (l*l + l2.*l2);
k = kern.variance*(2*(l*l2)./w2).^(kern.inputDimension/2).*exp(-n2./w2);
for i = 1:size(x, 2)
  gX(:, i) = 2*(X2(:, i) - x(i))./w2.*k;
end

base2 = k.*(kern.inputDimension/2*(l2.^4 - l^4) + 2*l*l*n2)./(w2.*w2.*l);
for i = 1:size(x, 2)
  gX(:, i) = gX(:, i) + base2*gl(i);
end
%base2 = k.*(kern.inputDimension/2*(l^4-l2.^4) + 2*l2.*l2.*n2)./ ...
%        (w2.*w2.*l2);

%for i = 1:size(x, 2)
%  gX(:, i) = gX(:, i) + sum(base2.*gl2(:, i));
%end
