function model = gpComputeAlpha(model, m)

% GPCOMPUTEALPHA Update the vector `alpha' for computing posterior mean quickly.
% FORMAT
% DESC updates the vectors that are known as `alpha' in the support
% vector machine, in other words invK*y, where y is the target values.
% ARG model : the model for which the alphas are going to be
% updated.
% ARG m : the values of m for which the updates will be made.
% RETURN model : the model with the updated alphas.
%
% SEEALSO : gpCreate, gpUpdateAD, gpUpdateKernels
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2009

% GP

if nargin < 2
  m = model.m;
end

switch model.approx
 case 'ftc'
  model.alpha = zeros(model.N, model.d);
  if ~isfield(model, 'isSpherical') | model.isSpherical
    model.alpha = model.invK_uu*m;
  else
    for i = 1:model.d
      ind = gpDataIndices(model, i);
      model.alpha(ind, i) = model.invK_uu{i}* ...
          m(ind, i);
    end
  end
 
 case {'dtc', 'dtcvar'}
  model.alpha = zeros(model.k, model.d);
  if ~isfield(model, 'isSpherical') | model.isSpherical
    model.alpha = model.Ainv*model.K_uf*m;
  else
    for i = 1:model.d
      ind = gpDataIndices(model, i);
      model.alpha(:, i) = model.Ainv{i} ...
          *model.K_uf(:, ind) ...
          *m(ind, i);
    end
  end
 case 'fitc'
  model.alpha = zeros(model.k, model.d);
  if ~isfield(model, 'isSpherical') | model.isSpherical
    model.alpha = model.Ainv*model.K_uf*model.Dinv*m;
  else
    for i = 1:model.d
      ind = gpDataIndices(model, i);
      model.alpha(:, i) = model.Ainv{i} ...
          *model.K_uf(:, ind) ...
          *model.Dinv{i}*m(ind, i);
    end
  end
 case 'pitc'
  model.alpha = zeros(model.k, model.d);
  if ~isfield(model, 'isSpherical') | model.isSpherical
    for i = 1:length(model.blockEnd)
      ind = gpBlockIndices(model, i);
      model.alpha = model.alpha+model.Ainv*model.K_uf(:, ind)* ...
          model.Dinv{i}*m(ind, :);
    end
  else
    for i = 1:length(model.blockEnd)
      for j = 1:model.d
        ind = gpDataIndices(model, j, i);
        model.alpha(:, j) = model.alpha(:, j)+model.Ainv{j}*model.K_uf(:, ind)* ...
            model.Dinv{i, j}*m(ind, j);
      end
    end  
  end
end
