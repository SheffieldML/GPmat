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

error('ouKernGradX not yet implemented.')