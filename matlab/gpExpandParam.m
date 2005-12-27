function model = gpExpandParam(model, params)

% GPEXPANDPARAMS Expand a parameter vector into a GP model.

% FGPLVM

startVal = 1;
endVal = model.k*model.q;
model.X_u = reshape(params(startVal:endVal), model.k, model.q);
startVal = endVal +1;
endVal = endVal + model.kern.nParams;
model.kern = kernExpandParam(model.kern, params(startVal:endVal));

switch model.approx
 case 'ftc'
  model = gpUpdateKernels(model, model.X, model.X_u);
 case {'dtc', 'fitc', 'pitc'}
  model = gpUpdateKernels(model, model.X, model.X_u, params(end));
 otherwise
  error('Unknown approximation type.')
end

if isfield(model, 'alpha')
  model = gpComputeAlpha(model);
end



