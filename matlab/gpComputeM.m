function m = gpComputeM(model)

% GPCOMPUTEM Compute the matrix m given the model.
%
% m = gpComputeM(model)
%

% Copyright (c) 2006 Neil D. Lawrence
% gpComputeM.m version 1.1



m = zeros(size(model.y));
for i = 1:model.d
  m(:, i) = (model.y(:, i) - model.bias(i));
  if model.scale(i)
    m(:, i) = m(:, i)/model.scale(i);
  end
end
