function obj = fmvuObjective(params, model)

% FMVUOBJECTIVE Log likelihood of FMVU model.
% FORMAT
% DESC computes the log likelihood of  the fast maximum variance unfolding model.
% ARG params : the parameters at which the model objective is evaluated.
% ARG model : the model structure for which the objective is being computed.
% RETURN obj : the computed objective function.
%
% SEEALSO : fmvuCreate, fmvuLogLikeGradients, modelObjective
%
% COPYRIGHT : Neil D. Lawrence, 2009

% MLTOOLS

  model = fmvuExpandParam(model, params);
  obj = sum(sum(model.kappa.*abs(model.delta2 - model.D2)));
  %obj = obj + sum(sum(model.X.*model.X));
end