function model = gpComputeAlpha(model, m)

% GPCOMPUTEALPHA Update the vector `alpha' for computing posterior mean quickly.
%
% model = gpComputeAlpha(model, m)
%

% Copyright (c) 2006 Neil D. Lawrence
% gpComputeAlpha.m version 1.4


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
 
 case 'dtc'
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

 case 'nftc'
  model.alpha = zeros(model.N, model.d);
  if ~isfield(model, 'isSpherical') | model.isSpherical
    model.alpha = model.Ainv*m;
  else
    for i = 1:model.d
      ind = gpDataIndices(model, i);
      model.alpha(ind, i) = model.Ainv{i}* ...
          m(ind, i);
    end
  end
end