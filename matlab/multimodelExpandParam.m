function model = multimodelExpandParam(model, params)

% MULTIMODELEXPANDPARAM Create model structure from MULTIMODEL model's parameters.
% FORMAT
% DESC returns a multi-task learning wrapper model structure filled with the
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
% SEEALSO : multimodelCreate, multimodelExtractParam, modelExpandParam
%
% COPYRIGHT : Neil D. Lawrence, 2007, 2008

% MLTOOLS
  
  endVal = model.numParams - model.numModels*model.numSep;
  baseParams = params(1:endVal);
  passParams(model.sharedIndices) = baseParams;
  if ~isempty(model.separateIndices)
    for i = 1:length(model.comp)
      startVal = endVal + 1
      endVal = endVal + model.numSep;
      passParams(model.separateIndices) = params(startVal:endVal);
      model.comp{i} = modelExpandParam(model.comp{i}, passParams);
    end
  end
end