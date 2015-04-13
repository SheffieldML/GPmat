function g = fmvuGradient(params, model)

% FMVUGRADIENT Gradient of FMVU model log likelihood with respect to parameters.
% FORMAT
% DESC computes the gradient of the fast maximum variance unfolding
% model's log likelihood with respect to the parameters.
% ARG model : model structure for which gradients are being
% computed.
% RETURN g : the returned gradients. 
%
% SEEALSO fmvuCreate, fmvuLogLikelihood, modelGradient 
%
% COPYRIGHT : Neil D. Lawrence, 2009

% MLTOOLS

  model = fmvuExpandParam(model, params);
  
  gX = model.L*model.X; % + model.X;
  gLambda = model.kappa.*abs(model.delta2 - model.D2);
  g = [2*gX(:)' gLambda(:)'];
end
