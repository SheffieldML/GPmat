function [model, ll] = gpsimMapUpdateKernels(model)

% GPSIMMAPUPDATEKERNELS Updates the kernel representations in the GPSIMMAP structure.
% FORMAT
% DESC updates any representations of the kernel in the model
% structure, such as invK, logDetK or K.
% ARG model : the model structure for which kernels are being
% updated.
% RETURN model : the model structure with the kernels updated.
% RETURN ll : the log likelihood of the model after update.
%
% SEEALSO gpsimMapExpandParam, gpsimMapCreate
%
% COPYRIGHT Neil D. Lawrence, 2006

% GPSIM

model.K = kernCompute(model.kern, model.mapt);
% Add jitter.
[model.invK, U, jitter] = pdinv(model.K);
if jitter>1e-4
  fprintf('Warning: gpsimMapUpdateKernels added jitter of %2.4f\n', jitter)
end
model.logDetK = logdet(model.K, U);
model = gpsimMapUpdatePosteriorCovariance(model);

if model.updateF
  options = defaultOptions;
  [model, ll] = gpsimMapUpdateF(model, options);
end