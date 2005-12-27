function model = gpComputeAlpha(model)

% GPCOMPUTEALPHA Update the vector `alpha' for computing posterior mean quickly.

% FGPLVM

switch model.approx
 case 'ftc'
  model.alpha = model.invK_uu*model.Y;
 
 case 'dtc'
  model.alpha = model.Ainv*model.K_uf*model.Y;
  
 case 'fitc'
  model.alpha = model.Ainv*model.K_uf*model.Dinv*model.Y;
  
 case 'pitc'
  model.alpha = zeros(size(model.X_u, 1), size(model.Y, 2));
  startVal = 1;
  for i = 1:length(model.blockEnd)
    endVal = model.blockEnd(i);
    ind = startVal:endVal;
    model.alpha = model.alpha+model.Ainv*model.K_uf(:, ind)* ...
        model.Dinv{i}*model.Y(ind, :);
    startVal = endVal + 1;
  end
  
end