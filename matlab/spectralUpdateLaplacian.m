function model = spectralUpdateLaplacian(model)
  
% SPECTRALUPDATELAPLACIAN Update the Laplacian using graph connections.
% FORMAT
% DESC updates the Laplacian matrix using values of the internode
% connections.
% ARG model : the model to be updated.
% RETURN model : the updated model.
% 
% COPYRIGHT : Neil D. Lawrence, 2009
%
% SEEALSO : spectralUpdateX, leOptimise

% MLTOOLS

for i = 1:size(model.indices, 1)
  for j = 1:model.k
    model.L(i, model.indices(i, j)) = -model.kappa(i, j);
    model.L(model.indices(i, j), i) = model.L(i, model.indices(i, j));
  end
end
% Set diagonal
D = sum(model.kappa, 2);
model.L(1:model.N+1:end) = D;
if model.isNormalised
  sqrtD = sqrt(D);
  model.L = (model.L.*(sqrtD*sqrtD'));
end
