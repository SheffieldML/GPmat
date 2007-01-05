function [params, names] = mlpExtractParam(model)

% MLPEXTRACTPARAM Extract weights and biases from an MLP.
% FORMAT
% DESC returns a vector of all the weights and biases from a
% multi-layer perceptron model. 
% ARG model : the model from which we wish to extract the weights
% and biases.
% RETURN params : vector of all the weights and biases returned by
% the model. The structure is governed by mlppak.
% RETURN names : optional additional returned cell array of the
% names of the parameters.
%
% SEEALSO : mlppak, mlpCreate, mlpExpandParam, modelExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2006

% MLTOOLS

params = mlppak(model);

if nargout > 1
  counter = 0;
  for j = 1:size(model.w1, 2)
    for i = 1:size(model.w1, 1)
      counter = counter + 1;
      names{counter} = ['Input weight ' num2str(i) '-' num2str(j)];
    end
  end
  for j = 1:size(model.b1, 2)
    counter = counter + 1;
    names{counter} = ['Hidden node bias ' num2str(j)];
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