function model = gpsimMapExpandParam(model, params)

% GPSIMMAPEXPANDPARAM Expand the given parameters into a GPSIMMAP structure.
% FORMAT
% DESC takes the given vector of parameters and places them in the
% model structure, it then updates any stored representations that
% are dependent on those parameters, for example kernel matrices
% etc..
% ARG model : the model structure for which parameters are to be
% updated.
% ARG params : a vector of parameters for placing in the model
% structure.
% RETURN model : a returned model structure containing the updated
% parameters.
% 
% SEEALSO : gpsimMapCreate, gpsimMapExtractParam, modelExtractParam, gpsimMapUpdateKernels
%
% COPYRIGHT : Neil D. Lawrence, 2006

% GPSIM
params = real(params);
if isfield(model, 'fix')
  for i = 1:length(model.fix)
    params(model.fix(i).index) = model.fix(i).value;
  end
end

if length(params) ~= model.numParams
  error('Parameter vector is incorrect length');
end
startVal = 1;
endVal = model.kern.nParams;
model.kern = kernExpandParam(model.kern, params(startVal:endVal));

% Check if there is a mean function.
if isfield(model, 'meanFunction') & ~isempty(model.meanFunction)
  startVal = endVal + 1;
  endVal = endVal + model.meanFunction.numParams;
  model.meanFunction = modelExpandParam(model.meanFunction, ...
                                        params(startVal:endVal));
end


model = gpsimMapUpdateKernels(model);
model = gpsimMapUpdatePosteriorCovariance(model);
model.updateW = true; % Needs to be re-done after parameter change.