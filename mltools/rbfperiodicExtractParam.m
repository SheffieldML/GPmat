function [params, names] = rbfperiodicExtractParam(model)

% RBFPERIODICEXTRACTPARAM Extract parameters from the RBFPERIODIC model structure.
% FORMAT
% DESC extracts parameters from the periodic radial basis function
% model structure into a vector of parameters for optimisation.
% ARG model : the model structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the model. 
%
% FORMAT
% DESC extracts parameters and parameter names from the periodic radial basis function
% model structure.
% ARG model : the model structure containing the parameters to be
% extracted.
% RETURN param : vector of parameters extracted from the model. 
% RETURN names : cell array of strings containing names for each
% parameter.
%
% SEEALSO rbfperiodicCreate, rbfperiodicExpandParam, modelExtractParam, scg, conjgrad
%
% COPYRIGHT : Neil D. Lawrence, 2007
%
% MLTOOLS

fhandle = str2func([model.widthTransform.type 'Transform']);
params = [model.thetaBar(:)' ...
          fhandle(model.sigma2(:)', 'xtoa') ...
          model.weights(:)'  ...
          model.bias(:)'];
if nargout > 1
  counter = 0;
  for j = 1:size(model.thetaBar, 2)
    for i = 1:size(model.thetaBar, 1)
      counter = counter + 1;
      names{counter} = ['Center position ' num2str(i) '-' num2str(j)];
    end
  end
  for j = 1:size(model.sigma2, 2)
    counter = counter + 1;
    names{counter} = ['Basis function width ' num2str(j)];
  end
  for j = 1:size(model.weights, 2)
    for i = 1:size(model.weights, 1)
      counter = counter + 1;
      names{counter} = ['Output weight ' num2str(i) '-' num2str(j)];
    end
  end
  for j = 1:size(model.bias, 2)
    counter = counter + 1;
    names{counter} = ['Output node bias ' num2str(j)];
  end
end
