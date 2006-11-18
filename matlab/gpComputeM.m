function m = gpComputeM(model)

% GPCOMPUTEM Compute the matrix m given the model.
% FORMAT
% DESC computes the matrix m (the scaled and bias removed matrix of
% the targets), given the model.
% ARG model : the model for which the values are to be computed.
% ARG m : the scaled and bias removed values.
%
% SEEALSO : gpCreate, gpComputeAlpha, gpUpdateAD
%
% COPYRIGHT : Neil D. Lawrence, 2006

% GP 

m = zeros(size(model.y));
for i = 1:model.d
  m(:, i) = (model.y(:, i) - model.bias(i));
  if model.scale(i)
    m(:, i) = m(:, i)/model.scale(i);
  end
end
