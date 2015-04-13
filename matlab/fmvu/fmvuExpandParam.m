function model = fmvuExpandParam(model, params)

% FMVUEXPANDPARAM Create model structure from FMVU model's parameters.
% FORMAT
% DESC returns a fast maximum variance unfolding model structure filled with the
% parameters in the given vector. This is used as a helper function to
% enable parameters to be optimised in, for example, the NETLAB
% optimisation functions.
% ARG model : the model structure in which the parameters are to be
% placed.
% ARG param : vector of parameters which are to be placed in the
% model structure.
% RETURN model : model structure with the given parameters in the
% relevant locations.
%
% SEEALSO : fmvuCreate, fmvuExtractParam, modelExpandParam
%
% COPYRIGHT : Neil D. Lawrence, 2009

% MLTOOLS

  startInd = 1;
  endInd = model.q*model.N;
  model.X = reshape(params(startInd:endInd), model.N, model.q);
  % Update the latent distances.
  for i = 1:size(model.indices, 1)
    for j = 1:model.k
      model.delta2(i, j) = dist2(model.X(i, :), model.X(model.indices(i, j), :));
    end
  end
  startInd = endInd + 1;
  endInd = endInd + model.k*model.N;
  model.kappa = reshape(exp(params(startInd:endInd)), model.N, model.k);
  signDiff = sign(model.delta2 - model.D2);
  model.kappa = model.kappa.* signDiff;
  % Update the stiffness matrix/Laplacian.
  model = spectralUpdateLaplacian(model);
  model.kappa = model.kappa.* signDiff;

end  
