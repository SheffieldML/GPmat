function model = pmvuExpandParam(model, params)

% PMVUEXPANDPARAM Create model structure from PMVU model's parameters.
% FORMAT
% DESC returns a probabilistic maximum variance unfolding model structure filled with the
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
% SEEALSO : pmvuCreate, pmvuExtractParam, modelExpandParam
%
% COPYRIGHT : Neil D. Lawrence 2009

% MLTOOLS

  fhandle = str2func([model.kappaTransform 'Transform']);
  params = fhandle(params, 'atox');
  model.kappa = reshape(params, model.N, model.k);
  model = spectralUpdateLaplacian(model);
  Kinv = model.L + model.gamma*speye(model.N);
  [model.K, U] = pdinv(full(Kinv));
  model.logDetK = - logdet(Kinv, U);
end
