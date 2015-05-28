function model = lfmUpdateKernels(model)

% LFMUPDATEKERNELS Updates the kernel representations in the LFM structure.
% FORMAT
% DESC updates any representations of the kernel in the model
% structure, such as invK, logDetK or K.
% ARG model : the model structure for which kernels are being
% updated.
% RETURN model : the model structure with the kernels updated.
%
% SEEALSO lfmExpandParam, lfmCreate
%
% COPYRIGHT Neil D. Lawrence, 2007

% KERN


model.K = real(kernCompute(model.kern, model.t)); 
if isfield(model, 'yvar')
  model.K = model.K + diag(model.yvar);
end
[model.invK, U, jitter] = pdinv(model.K);
if jitter>1e-4
  fprintf('Warning: lfmUpdateKernels added jitter of %2.4f\n', jitter)
end
model.logDetK = logdet(model.K, U);
