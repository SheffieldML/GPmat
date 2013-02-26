function model = ivmUpdateM(model, index)

% IVMUPDATEM Update matrix M, L, v and mu.
%
%	Description:
%
%	MODEL = IVMUPDATEM(MODEL, INDEX) updates the stored representations
%	in the IVM (M, L, v and mu) given a new data point.
%	 Returns:
%	  MODEL - the returned model with all the representations up to
%	   date.
%	 Arguments:
%	  MODEL - the model for which the represenations are to be updated.
%	  INDEX - the index of the data point that is being included.
%	
%
%	See also
%	IVMCREATE, IVMADDPOINT, KERNCOMPUTE


%	Copyright (c) 2005 Neil D. Lawrence

activePoint = length(model.I)+1;

% Compute the kernel(s) at the new point.
model.kern.Kstore(:, activePoint) = kernCompute(model.kern, model.X, ...
                                                    model.X(index, :));
if isfield(model.kern, 'whiteVariance')
  model.kern.Kstore(index, activePoint) = model.kern.Kstore(index, activePoint) ...
      + model.kern.whiteVariance;
end
% Compute the values of M for the new point
for c = 1:length(model.Sigma)
  if length(model.I) > 0
    a = model.Sigma(c).M(:, index);
    s = model.kern.Kstore(:, activePoint)' - a'*model.Sigma(c).M;
    lValInv = sqrt(model.nu(index, c));
    model.Sigma(c).M = [model.Sigma(c).M; lValInv*s];
    ainv = (-a*lValInv)'/model.Sigma(c).L;
    model.Sigma(c).L = [[model.Sigma(c).L; a'] [zeros(length(model.I),1); 1/lValInv]];
    model.Sigma(c).Linv = [[model.Sigma(c).Linv; ainv] [zeros(length(model.I),1); lValInv]];
  else
    s = model.kern.Kstore(:, 1)';
    lValInv = sqrt(model.nu(index, c));
    model.Sigma(c).M = [lValInv*s];
    
    model.Sigma(c).L = 1/lValInv;
    model.Sigma(c).Linv = lValInv;
  end
    
  model.varSigma(:, c) = model.varSigma(:, c) - ((model.nu(index, c)*s).*s)';
  model.mu(:, c) = model.mu(:, c) + model.g(index, c)*s';  
end
if length(model.Sigma)==1 & size(model.y, 2)>1
  for c = 2:size(model.y, 2)
    model.varSigma(:, c) = model.varSigma(:, c) ...
        - ((model.nu(index, c)*s).*s)';
    model.mu(:, c) = model.mu(:, c) + model.g(index, c)*s'; 
  end
end
