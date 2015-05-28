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

  model.L = spalloc(model.N, model.N, model.N*model.k);
  for i = 1:size(model.indices, 1)
    for j = 1:model.k
      k = model.indices(i, j);
      if i>k
        model.L(i, k) = model.L(i, k)-model.kappa(i, j);        
        model.L(k, i) = model.L(i, k); 
      else
        model.L(k, i) = model.L(k, i)-model.kappa(i, j);
        model.L(i, k) = model.L(k, i); 
      end
    end
  end
  
  % Set diagonal
  D = sum(model.kappa, 2);
  D = -sum(model.L, 2);
  model.L(1:model.N+1:end) = D;
  if model.isNormalised
    sqrtD = sqrt(D);
    model.L = (model.L.*((1./sqrtD)*(1./sqrtD')));
  end
end
