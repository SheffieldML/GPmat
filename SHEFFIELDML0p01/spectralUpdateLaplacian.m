function model = spectralUpdateLaplacian(model)

% SPECTRALUPDATELAPLACIAN Update the Laplacian using graph connections.
%
%	Description:
%
%	MODEL = SPECTRALUPDATELAPLACIAN(MODEL) updates the Laplacian matrix
%	using values of the internode connections.
%	 Returns:
%	  MODEL - the updated model.
%	 Arguments:
%	  MODEL - the model to be updated.
%	
%
%	See also
%	SPECTRALUPDATEX, LEOPTIMISE


%	Copyright (c) 2009 Neil D. Lawrence


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
