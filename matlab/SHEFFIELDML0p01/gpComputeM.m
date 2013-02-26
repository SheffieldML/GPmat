function m = gpComputeM(model)

% GPCOMPUTEM Compute the matrix m given the model.
%
%	Description:
%
%	M = GPCOMPUTEM(MODEL) computes the matrix m (the scaled, bias and
%	mean function removed matrix of the targets), given the model.
%	 Returns:
%	  M - the scaled, bias and mean function removed values.
%	 Arguments:
%	  MODEL - the model for which the values are to be computed.
%	
%
%	See also
%	GPCREATE, GPCOMPUTEALPHA, GPUPDATEAD


%	Copyright (c) 2006 Neil D. Lawrence


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
