function gX = diagKernGradX(kern, X, X2)

% DIAGKERNGRADX Gradient of DIAG kernel with respect to a point x.
% FORMAT
% DESC computes the gradient of the diagonal noise covariance function
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
% DESC computes the gradident of the diagonal noise covariance function
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
% SEEALSO diagKernParamInit, kernGradX, diagKernDiagGradX
%
% COPYRIGHT : Neil D. Lawrence, 2011

% KERN

  if size(X, 2)>1
    error('Diag kernel requires 1-dimensional input.')
  end
  gX = zeros(size(X2, 1), 1, size(X, 1));
  trans = str2func([kern.trans, 'Transform']);
  vars = trans(X, 'atox');
  gXtemp = kern.variance*trans(vars, 'gradfact');
  for i = 1:size(X, 1);
    gX(i, :, i) = 0.5*gXtemp(i);
  end
end
