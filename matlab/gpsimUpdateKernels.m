function model = gpsimUpdateKernels(model)

% GPSIMUPDATEKERNELS Updates the kernel representations in the GPSIM structure.
% FORMAT
% DESC updates any representations of the kernel in the model
% structure, such as invK, logDetK or K.
% ARG model : the model structure for which kernels are being
% updated.
% RETURN model : the model structure with the kernels updated.
%
% SEEALSO gpsimExpandParam, gpsimCreate
%
% COPYRIGHT Neil D. Lawrence, 2006
%  
% MODIFIED : Pei Gao, 2008

% SHEFFIELDML
  
eps = 1e-6;                             % minimum noise variance for the
                                        % RBF kernel of TF.

if isfield(model, 'proteinPrior') && isfield(model, 'timesCell')
  k = real(kernCompute(model.kern, model.timesCell));
  if model.includeNoise
    noiseVar = [zeros(model.kern.comp{1}.diagBlockDim{1}, 1); model.yvar];
  else
    noiseVar = [eps*ones(model.kern.diagBlockDim{1}, 1); model.yvar];
  end
else
  k = real(kernCompute(model.kern, model.t));
  noiseVar = model.yvar;
end

model.K = k + diag(noiseVar);
[model.invK, U, jitter] = pdinv(model.K);
if jitter>1e-4
  warning('gpsim:jitter', 'gpsimUpdateKernels added jitter of %2.4f\n', jitter)
end
model.logDetK = logdet(model.K, U);
