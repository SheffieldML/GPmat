function model = ivmReconstruct(kern, noise, ivmInfo, X, y)

% IVMRECONSTRUCT Reconstruct an IVM form component parts.
% FORMAT
% DESC takes component parts of an IVM model and reconstructs the
% IVM model. The component parts are normally retrieved from a
% saved file.
% ARG kern : a kernel structure for the IVM.
% ARG noise : a noise structure for the IVM.
% ARG ivmInfo : the active set and the inactive set of the IVM as
% well as the site parameters, stored in a structure.
% ARG X : the input training data for the IVM.
% ARG y : the output target training data for the IVM.
% RETURN model : an IVM model structure that combines the component
% parts.
% 
% SEEALSO : ivmDeconstruct, ivmCreate, ivmComputeLandM
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

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
