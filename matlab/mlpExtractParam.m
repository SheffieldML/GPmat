function [params, names] = mlpExtractParam(model)

% MLPEXTRACTPARAM Extract weights and biases from an MLP.
% FORMAT
% DESC returns a vector of all the weights and biases from a
% multi-layer perceptron model. For single hidden layer models the
% function is a wrapper for the mlppak command.
% ARG model : the model from which we wish to extract the weights
% and biases.
% RETURN params : vector of all the weights and biases returned by
% the model. The structure is governed by mlppak.
% RETURN names : optional additional returned cell array of the
% names of the parameters.
%
% SEEALSO : mlppak, mlpCreate, mlpExpandParam, modelExtractParam
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2007

% MLTOOLS

if length(model.hiddenDim) == 1
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
else
  params = zeros(1, model.numParams);
  startVal = 1;
  endVal = model.inputDim*model.hiddenDim(1);
  params(startVal:endVal) = model.w{1}(:)';
  startVal = endVal + 1;
  endVal = endVal + model.hiddenDim(1);
  params(startVal:endVal) = model.b{1};
  for i = 2:length(model.hiddenDim)
    startVal = endVal + 1;
    endVal = endVal + model.hiddenDim(i-1)*model.hiddenDim(i);
    params(startVal:endVal) = model.w{i}(:)';
    startVal = endVal + 1;
    endVal = endVal + model.hiddenDim(i);
    params(startVal:endVal) = model.b{i};
  end
  i = length(model.hiddenDim);
  startVal = endVal + 1;
  endVal = endVal + model.hiddenDim(i)*model.outputDim;
  params(startVal:endVal) = model.w{i+1}(:)';
  startVal = endVal + 1;
  endVal = endVal + model.outputDim;
  params(startVal:endVal) = model.b{i+1};
  if nargout > 1
    counter = 0;
    for j = 1:size(model.w{1}, 2)
      for i = 1:size(model.w{1}, 1)
        counter = counter + 1;
        names{counter} = ['Input weight ' num2str(i) '-' num2str(j)];
      end
    end
    for j = 1:size(model.b{1}, 2)
      counter = counter + 1;
      names{counter} = ['Hidden node bias ' num2str(j)];
    end
    for k = 2:length(model.hiddenDim)
      for j = 1:size(model.w{k}, 2)
        for i = 1:size(model.w{k}, 1)
          counter = counter + 1;
          names{counter} = ['Hidden weight layer ' num2str(k-1) '-' ...
                            num2str(k) ', node '  num2str(i) '-' ...
                            num2str(j)];
        end
      end
    end
    for j = 1:size(model.w{end}, 2)
      for i = 1:size(model.w{end}, 1)
        counter = counter + 1;
        names{counter} = ['Output weight ' num2str(i) '-' num2str(j)];
      end
    end
    for j = 1:size(model.b2, 2)
      counter = counter + 1;
      names{counter} = ['Output node bias ' num2str(j)];
    end
  end
end
