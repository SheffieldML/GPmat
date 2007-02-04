function g = multimodelLogLikeGradients(model)

% MULTIMODELLOGLIKEGRADIENTS Gradient of MULTIMODEL model log likelihood with respect to parameters.
% FORMAT
% DESC computes the gradient of the multi-task learning wrapper
% model's log likelihood with respect to the parameters.
% ARG model : model structure for which gradients are being
% computed.
% RETURN g : the returned gradients. 
%
% SEEALSO multimodelCreate, multimodelLogLikelihood, modelLogLikeGradients 
%
% COPYRIGHT : Neil D. Lawrence, 2007

% MLTOOLS

g = zeros(1, model.numParams);
for i = 1:length(model.comp)
  g = g + modelLogLikeGradients(model.comp{i});
end