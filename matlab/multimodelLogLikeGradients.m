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
%
% MODIFICATIONS : Mauricio Alvarez, 2009

% MLTOOLS

  g = zeros(1, model.numParams);
  if numel(model.outputDim) == 1
  endShared = model.numParams - model.numModels*model.numSep;
  endVal = endShared; 
  for i = 1:length(model.comp)
    gModel = modelLogLikeGradients(model.comp{i}); 
    g(1:endShared) = g(1:endShared) + gModel(model.sharedIndices);
    if ~isempty(model.separateIndices)
      startVal = endVal + 1;
      endVal = endVal + model.numSep;
      %g(startVal:endVal) = g(model.separateIndices);
      g(startVal:endVal) = gModel(model.separateIndices);      
    end
  end
  else
      startVal = 1;
      endVal = 0;
      for i = 1:length(model.comp)
          endVal = endVal + model.comp{i}.nParams;
          if (i==6) || (i ==14)
              a =1;
          end
          g(startVal:endVal) = modelLogLikeGradients(model.comp{i});
          startVal = endVal + 1;
      end
  end
end
