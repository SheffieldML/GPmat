function model = gpUpdateKernels(model, X, X_u, rawBeta)

% GPUPDATEKERNELS Update the kernels that are needed.
% FORMAT
% DESC updates any representations of the kernel in the model
% structure, such as invK, logDetK or K.
% ARG model : the model structure for which kernels are being
% updated.
% RETURN model : the model structure with the kernels updated.
%
% SEEALSO gpExpandParam, gpCreate
%
% COPYRIGHT Neil D. Lawrence, 2005, 2006

% GP

jitter = 1e-6;

if nargin < 4
  rawBeta = [];
end
if ~isempty(rawBeta)
  fhandle = str2func([model.betaTransform 'Transform']);
  model.beta = fhandle(rawBeta, 'atox');
end
switch model.approx
 case 'ftc'
  % Long term should allow different kernels in each dimension here.
  model.K_uu = kernCompute(model.kern, X);
    
  if ~isfield(model, 'isSpherical') | model.isSpherical
    % Add inverse beta to diagonal if it exists.
    if isfield(model, 'beta')
      model.K_uu(1:size(model.K_uu, 1)+1:end) = ...
          model.K_uu(1:size(model.K_uu, 1)+1:end) + 1./model.beta;
    end
    [model.invK_uu, U] = pdinv(model.K_uu);
    model.logDetK_uu = logdet(model.K_uu, U);
  else   
    % TODO: Need to deal with adding in beta here.
    for i = 1:model.d
      ind = gpDataIndices(model, i);
      [model.invK_uu{i}, U] = pdinv(model.K_uu(ind, ind));
      model.logDetK_uu(i) = logdet(model.K_uu(ind, ind), U);
    end
  end
 case {'dtc', 'fitc', 'pitc'}
  model.K_uu = kernCompute(model.kern, X_u);
  
  if ~isfield(model.kern, 'whiteVariance') | model.kern.whiteVariance == 0
    % There is no white noise term so add some jitter.
    model.K_uu = model.K_uu ...
        + sparseDiag(repmat(jitter, size(model.K_uu, 1), 1));
  end
  model.K_uf = kernCompute(model.kern, X_u, X);
  [model.invK_uu, model.sqrtK_uu] = pdinv(model.K_uu);
  model.logDetK_uu = logdet(model.K_uu, model.sqrtK_uu);

end

model = gpUpdateAD(model);  
