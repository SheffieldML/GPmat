function model = gpReconstruct(kern, noise, gpInfo, X, y)

% GPRECONSTRUCT Reconstruct an GP form component parts.
% FORMAT
% DESC takes component parts of an GP model and reconstructs the
% GP model. The component parts are normally retrieved from a
% saved file.
% ARG kern : a kernel structure for the GP.
% ARG noise : a noise structure for the GP (currently ignored).
% ARG gpInfo : the active set and other information stored in a structure.
% ARG X : the input training data for the GP.
% ARG y : the output target training data for the GP.
% RETURN model : an GP model structure that combines the component
% parts.
% 
% SEEALSO : gpDeconstruct, gpCreate
%
% COPYRIGHT : Neil D. Lawrence, 2007, 2009

% GP

options = gpOptions(gpInfo.approx);
options.kern = kern;
switch gpInfo.approx
 case 'ftc'
 case {'dtc', 'dtcvar', 'fitc', 'pitc'}
  options.numActive = size(gpInfo.X_u, 1);
end
model = gpCreate(size(X, 2), size(y, 2), X, y, options);
model.scale = gpInfo.scale;
model.bias = gpInfo.bias;
model.m = gpComputeM(model);
model.learnScales = gpInfo.learnScales;
switch model.approx
 case 'ftc'
 case {'dtc', 'dtcvar', 'fitc', 'pitc'}
  model.beta = gpInfo.beta;
  model.fixInducing = gpInfo.fixInducing;
  if gpInfo.fixInducing
    model.inducingIndices = gpInfo.inducingIndices;
  else
    model.X_u = gpInfo.X_u;
  end
end

if gpInfo.d ~= size(y, 2)
  error('y does not have correct number of dimensions.')
end
if gpInfo.q ~= size(X, 2)
  error('X does not have correct number of dimensions.')
end
% FOrce update of everything.
params = gpExtractParam(model);
model = gpExpandParam(model, params);
