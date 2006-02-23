function model = ivmComputeLandM(model)

% IVMCOMPUTELANDM Compute the L and M matrix.

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
