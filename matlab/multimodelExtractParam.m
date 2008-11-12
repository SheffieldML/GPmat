function [passParams, passNames] = multimodelExtractParam(model)

% MULTIMODELEXTRACTPARAM Extract parameters from the MULTIMODEL model structure.
% FORMAT
% DESC extracts parameters from the multi-task learning wrapper
% model structure into a vector of parameters for optimisation.
% ARG model : the model structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the model. 
%
% FORMAT
% DESC extracts parameters and parameter names from the multi-task learning wrapper
% model structure.
% ARG model : the model structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the model. 
% RETURN names : cell array of strings containing names for each
% parameter.
%
% SEEALSO multimodelCreate, multimodelExpandParam, modelExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2007, 2008
%
% MLTOOLS

  passParams = zeros(1, model.numParams);
  passNames = cell(1, model.numParams);
  if nargout > 1
    [receiveParams, receiveNames] = modelExtractParam(model.comp{1});
  else
    receiveParams = modelExtractParam(model.comp{1});
  end
  endVal = model.numParams - model.numSep;
  passParams(1:endVal) = receiveParams(model.sharedIndices);
  if nargout > 1
    passNames{1:endVal} = receiveNames{model.sharedIndices};
  end
  if ~isempty(model.separateIndices)
    startVal = endVal + 1;
    endVal = endVal + model.numSep;
    passParams(startVal:endVal) = receiveParams(model.separateIndices);
    if nargout > 1
      passNames{startVal:endVal} = receiveNames{model.separateIndices};
    end
    for i = 2:length(model.comp)
      startVal = endVal+1
      endVal = endVal + model.numSep;
      if nargout > 1
        [receiveParams, receiveNames] = modelExtractParam(model.comp{2});
      else
        params = modelExtractParam(model.comp{2});
      end
      passParams(startVal:endVal) = receiveParams(model.separateIndices);
      if nargout > 1
        passNames{startVal:endVal} = receiveNames{model.separateIndices};
      end
    end
  end
end