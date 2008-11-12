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
% COPYRIGHT : Neil D. Lawrence, 2007, 2008

% MLTOOLS

  g = zeros(1, model.numParams);
  endShared = model.numParams - model.numModels*model.numSep;
  endVal = endShared; 
  for i = 1:length(model.comp)
    gModel = modelLogLikeGradients(model.comp{i}); 
    g(1:endShared) = g(1:endShared) + gModel(model.sharedIndices);
    if ~isempty(model.separateIndices)
      startVal = endVal + 1;
      endVal = endVal + model.numSep;
      g(startVal:endVal) = g(model.separateIndices);
    end
  end
end
