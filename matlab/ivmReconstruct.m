function model = ivmReconstruct(kern, noise, ivmInfo, X, y)

% IVMLOAD Load an IVM structure.

% IVM


model.X = X;
model.y = y;

model.kern = kern;
model.noise = noise;

model.I = ivmInfo.I;
model.d = length(model.I);
model.J = ivmInfo.J;

model.m = ivmInfo.m;
model.beta = ivmInfo.beta;

model.kern.Kstore = kernCompute(model.kern, model.X, ...
                                model.X(model.I, :));
if isfield(model.kern, 'whiteVariance')
  model.kern.Kstore(model.I, 1:model.d) = ...
      model.kern.Kstore(model.I, 1:model.d) ...
      + model.kern.whiteVariance;
end
model.kern.diagK = kernDiagCompute(model.kern, model.X);

if model.noise.spherical
  model.Sigma.L = chol(model.kern.Kstore(model.I, :) ...
                       + diag(1./model.beta(model.I)))';
  model.Sigma.Linv = eye(size(model.Sigma.L))/model.Sigma.L;
  model.Sigma.M = model.Sigma.Linv*model.kern.Kstore';
  model.Sigma.robust = 0;
else
  for i = 1:size(y, 2)
    model.Sigma(i).L = chol(model.kern.Kstore(model.I, :) ...
                            + diag(1./model.beta(model.I, i)))';
    model.Sigma(i).Linv = eye(size(model.Sigma(i).L))/model.Sigma(i).L;
    model.Sigma(i).M = model.Sigma(i).Linv*model.kern.Kstore';
    model.Sigma(i).robust = 0;
  end
end
model.d = length(model.I);
[model.mu, model.varSigma] = ivmPosteriorMeanVar(model, X);
