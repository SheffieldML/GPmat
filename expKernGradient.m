function g = expKernGradient(kern, varargin)

% EXPKERNGRADIENT Gradient of EXP kernel's parameters.
% FORMAT
% DESC computes the gradient of functions with respect to the
% exponentiated kernel's parameters. As well as the kernel
% structure and the input positions, the user provides a matrix
% PARTIAL which gives the partial derivatives of the function with
% respect to the relevant elements of the kernel matrix. 
% ARG kern : the kernel structure for which the gradients are being
% computed.
% ARG x : the input locations for which the gradients are being
% computed. 
% ARG partial : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The argument takes
% the form of a square matrix of dimension  numData, where numData is
% the number of rows in X.
% RETURN g : gradients of the function of interest with respect to
% the kernel parameters. The ordering of the vector should match
% that provided by the function kernExtractParam.
%
% FORMAT
% DESC computes the derivatives as above, but input locations are
% now provided in two matrices associated with rows and columns of
% the kernel matrix. 
% ARG kern : the kernel structure for which the gradients are being
% computed.
% ARG x1 : the input locations associated with the rows of the
% kernel matrix.
% ARG x2 : the input locations associated with the columns of the
% kernel matrix.
% ARG partial : matrix of partial derivatives of the function of
% interest with respect to the kernel matrix. The matrix should
% have the same number of rows as X1 and the same number of columns
% as X2 has rows.
% RETURN g : gradients of the function of interest with respect to
% the kernel parameters.
%
% SEEALSO expKernParamInit, kernGradient, expKernDiagGradient, kernGradX
%
% COPYRIGHT : Neil D. Lawrence, 2006

% KERN



[K, argK, kDiagMat] = expKernCompute(kern, varargin{1:end-1});
g = zeros(1, kern.nParams);
g(1) = sum(sum(K.*varargin{end}))/kern.variance;
if kern.isStationary
  % Modifiy the partial derivatives (to 'trick' kernGradient into
  % providing correct gradients)
  covGrad = kern.variance*varargin{end}.*exp(argK)*exp(kDiagMat);
  g(2:end) = kernGradient(kern.argument, varargin{1:end-1}, covGrad);
  % Modifiy the partial derivatives 
  covGrad = kern.variance*varargin{end}.*(exp(argK)-1)*exp(kDiagMat);
  g(2:end) = g(2:end) + kernDiagGradient(kern.argument, varargin{1}, covGrad(:));
else
  % Modifiy the partial derivatives (to 'trick' kernGradient into
  % providing correct gradients)
%   covGrad = kern.variance*varargin{end}.*exp(argK)*exp(kDiagMat);
%   g(2:end) = kernGradient(kern.argument, varargin{1:end-1}, covGrad);
  % Modifiy the partial derivatives 
%  covGrad = kern.variance*varargin{end}.*(exp(argK)-1)*exp(kDiagMat);
%  g(2:end) = g(2:end) + kernDiagGradient(kern.argument, varargin{1}, covGrad(:));
end
