function gX = diagKernDiagGradX(kern, X)

% DIAGKERNDIAGGRADX Gradient of DIAG kernel's diagonal with respect to X.
% FORMAT
% DESC computes the gradient of the diagonal of the diagonal noise covariance function kernel matrix with
% respect to the elements of the design matrix given in X.
% ARG kern : the kernel structure for which gradients are being computed.
% ARG X : the input data in the form of a design matrix.
% RETURN gX : the gradients of the diagonal with respect to each element
% of X. The returned matrix has the same dimensions as X.
%
% SEEALSO : diagKernParamInit, kernDiagGradX, diagkernGradX
%
% COPYRIGHT : Neil D. Lawrence, 2011

% KERN

  if size(X, 2)>1
    error('Diag kernel requires 1-dimensional input.')
  end
  trans = str2func([kern.trans, 'Transform']);
  vars = trans(X, 'atox');
  gX = kern.variance*trans(vars, 'gradfact');
end
