function m = gpComputeM(model)

% GPCOMPUTEM Compute the matrix m given the model.
% FORMAT
% DESC computes the matrix m (the scaled, bias and mean function
% removed matrix of the targets), given the model.
% ARG model : the model for which the values are to be computed.
% RETURN m : the scaled, bias and mean function removed values.
%
% SEEALSO : gpCreate, gpComputeAlpha, gpUpdateAD
%
% COPYRIGHT : Neil D. Lawrence, 2006

% GP 

% Remove mean function value from m (if mean function present).
if isfield(model, 'meanFunction') && ~isempty(model.meanFunction)
  m = model.y - modelOut(model.meanFunction, model.X);
else
  m = model.y;
end

% Remove bias and apply scale.
for i = 1:model.d
  m(:, i) = m(:, i) - model.bias(i);
  if model.scale(i)
    m(:, i) = m(:, i)/model.scale(i);
  end
end
