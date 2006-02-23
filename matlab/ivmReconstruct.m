function model = ivmReconstruct(kern, noise, ivmInfo, X, y)

% IVMRECONSTRUCT Reconstruct an IVM form component parts.

% IVM

model.type = 'ivm';

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

model = ivmComputeLandM(model);
model.d = length(model.I);
[model.mu, model.varSigma] = ivmPosteriorMeanVar(model, X);
