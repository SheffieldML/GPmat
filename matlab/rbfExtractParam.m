function [params, names] = rbfExtractParam(model)

% RBFEXTRACTPARAM Wrapper for NETLAB's rbfpak.
% FORMAT
% DESC returns a vector of all the weights and biases from a
% RBF network model. For single hidden layer models the
% function is a wrapper for the rbfpak command.
% ARG model : the model from which we wish to extract the weights
% and biases.
% RETURN params : vector of all the weights and biases returned by
% the model. The structure is governed by rbfpak.
% RETURN names : optional additional returned cell array of the
% names of the parameters.
%
% SEEALSO : rbfpak, rbfCreate, rbfExpandParam, modelExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2007, 2008


% MLTOOLS

params = rbfpak(model);
if nargout > 1
  counter = 0;
  for j = 1:size(model.c, 2)
    for i = 1:size(model.c, 1)
      counter = counter + 1;
      names{counter} = ['Input centre ' num2str(i) '-' num2str(j)];
    end
  end
  for j = 1:size(model.wi, 2)
    counter = counter + 1;
    names{counter} = ['Hidden node width ' num2str(j)];
  end
  for j = 1:size(model.w2, 2)
    for i = 1:size(model.w2, 1)
      counter = counter + 1;
      names{counter} = ['Output weight ' num2str(i) '-' num2str(j)];
    end
  end
  for j = 1:size(model.b2, 2)
    counter = counter + 1;
    names{counter} = ['Output node bias ' num2str(j)];
  end
end
