function model = ivmComputeLandM(model)

% IVMCOMPUTELANDM Compute the L and M matrix.
% FORMAT
% DESC computes the matrices L and M from a given active set and
% site parameters. Normally L and M are developed in a sequential
% way, as points are added. This function creates them directly
% from a given active set, for use, e.g. when loading in a saved
% model from file (ivmReconstruct).
% ARG model : the model for which we are recomputing L and M.
% RETURN model : the model with the updated L and M values.
%
% COPYRIGHT : Neil D. Lawrence, 2005
%
% SEEALSO : ivmReconstruct

% IVM

if model.noise.spherical
  model.Sigma.L = chol(model.kern.Kstore(model.I, :) ...
                       + diag(1./model.beta(model.I)))';
  model.Sigma.Linv = eye(size(model.Sigma.L))/model.Sigma.L;
  model.Sigma.M = model.Sigma.Linv*model.kern.Kstore';
else
  for i = 1:size(model.y, 2)
    model.Sigma(i).L = chol(model.kern.Kstore(model.I, :) ...
                            + diag(1./model.beta(model.I, i)))';
    model.Sigma(i).Linv = eye(size(model.Sigma(i).L))/model.Sigma(i).L;
    model.Sigma(i).M = model.Sigma(i).Linv*model.kern.Kstore';
  end
end
