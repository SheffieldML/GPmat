function model = gpComputeAlpha(model, m)

% GPCOMPUTEALPHA Update the vector `alpha' for computing posterior mean quickly.
% FGPLVM

if nargin < 2
  m = model.m;
end

switch model.approx
 case 'ftc'
  model.alpha = model.invK_uu*m;
 
 case 'dtc'
  model.alpha = model.Ainv*model.K_uf*m;
  
 case 'fitc'
  model.alpha = model.Ainv*model.K_uf*model.Dinv*m;
  
 case 'pitc'
  model.alpha = zeros(size(model.X_u, 1), size(m, 2));
  startVal = 1;
  for i = 1:length(model.blockEnd)
    endVal = model.blockEnd(i);
    ind = startVal:endVal;
    model.alpha = model.alpha+model.Ainv*model.K_uf(:, ind)* ...
        model.Dinv{i}*m(ind, :);
    startVal = endVal + 1;
  end
  
end