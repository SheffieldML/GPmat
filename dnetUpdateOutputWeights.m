function model = dnetUpdateOutputWeights(model)
  
% DNETUPDATEOUTPUTWEIGHTS Do an M-step (update parameters) on an Density Network model.
% FORMAT
% DESC updates the parameters of a Density Network model.
% ARG model : the model which is to be updated.
% RETURN model : the model with updated parameters.
%
% SEEALSO : dnetCreate, dnetUpdateBeta
%
% COPYRIGHT : Neil D. Lawrence, 2008

% MLTOOLS

  if model.basisStored
    G = sparseDiag(sum(model.w, 1)');
    Phi = [model.Phi ones(model.M, 1)];
    A = pdinv(Phi'*G*Phi+model.alpha/model.beta*speye(size(Phi, 2)))*Phi'*model.w'*model.y;
    model.A = A(1:size(model.Phi, 2), :);
    model.b = A(end, :);
    
    model.mapping = modelSetOutputWeights(model.mapping, model.A, model.b);
  end
