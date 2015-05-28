function gT = ouKernGradX(kern, t1, t2)

% OUKERNGRADX Gradient of OU kernel with respect to a point x (see
% ouKernCompute or ouKernParamInit for a more detailed description of the
% OU kernel).
% FORMAT
% DESC computes the gradient of the Ornstein-Uhlenbeck kernel
% kernel with respect to the input positions. 
% ARG kern : kernel structure for which gradients are being
% computed.
% ARG t1 : locations against which gradients are being computed.
% RETURN gT : the returned gradients. The gradients are returned in
% a matrix which is numData x numInputs x numData. Where numData is
% the number of data points and numInputs is the number of input
% dimensions in t1 (currently always one).
%
% FORMAT
% DESC computes the gradident of the Ornstein-Uhlenbeck kernel
% kernel with respect to the input positions where both the row
% positions and column positions are provided separately.
% ARG kern : kernel structure for which gradients are being
% computed.
% ARG t1 : row locations against which gradients are being computed.
% ARG t2 : column locations against which gradients are being computed.
% RETURN gT : the returned gradients. The gradients are returned in
% a matrix which is numData2 x numInputs x numData1. Where numData1 is
% the number of data points in t1, numData2 is the number of data
% points in t2 and numInputs is the number of input dimensions in t1
% (currently always one).
%
% SEEALSO ouKernParamInit, kernGradX, ouKernDiagGradX
%
% COPYRIGHT : David Luengo, 2009

% KERN

if nargin < 3
  t2 = t1;
end
if size(t1, 2) > 1 | size(t2, 2) > 1
  error('Input can only have one column');
end

gT = zeros(size(t1, 1), 1, size(t2, 1));

c = - 0.5 * kern.variance;
decay = kern.decay;
isStationary = kern.isStationary;
for i = size(t1, 1)
    % Gradient when t1 < t2 (i.e. \Delta_t < 0)
    s = sign(t1(i)-t2);
    gTemp = s .* exp(-decay*abs(t1(i)-t2));
    if (isStationary == false)
        gTemp =  gTemp - exp(-decay*(t1(i)+t2));
    end
    gT(i, 1, :) = c * gTemp;
end
