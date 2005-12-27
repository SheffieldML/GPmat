function model = fgplvmExpandParam(model, params)

% FGPLVMEXPANDPARAM Expand a parameter vector into a GP-LVM model.

% FGPLVM

startVal = 1;
if isfield(model, 'back')
  endVal = model.back.numParams;
  model.back = modelExpandParam(model.back, params(startVal:endVal));
  model.X = modelOut(model.back, model.Y);
else
  endVal = model.N*model.q;
  model.X = reshape(params(startVal:endVal), model.N, model.q);
end
startVal = endVal+1;
endVal = endVal + model.k*model.q + model.kern.nParams;

switch model.approx
 case 'ftc'
  endVal = endVal;
 case {'dtc', 'fitc', 'pitc'}
  endVal = endVal + 1; % account for sigma2 attached to the end.
 otherwise
  error('Unknown approximation type.')
end

model = gpExpandParam(model, params(startVal:endVal));

if isfield(model, 'dynamics') & ~isempty(model.dynamics)
  model.dynamics = gpUpdateKernels(model.dynamics, model.X(1:end-1, :), ...
                                   model.X_u);    
end


