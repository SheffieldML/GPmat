function g = dnetOutputGrad(model, X)

% DNETOUTPUTGRAD Evaluate derivatives of dnet model outputs with respect to parameters.
% FORMAT
% DESC evaluates the derivates of a density network's
% outputs with respect to the parameters of the mapping function.
% ARG model : the model for which the derivatives are to be
% computed.
% ARG X : the input data locations where the gradients are to be
% computed.
% RETURN g : the gradient of the outputs of the density network
% perceptron with respect to each of the parameters. The size of
% the matrix is number of data x number of parameters x number of
% outputs of the model.
%
% SEEALSO : modelOutputGrad, dnetCreate
%
% COPYRIGHT : Neil D. Lawrence, 2008

% MLTOOLS

  g = modelOutputGrad(model.mapping, X);
  % Add output gradient wrt beta.
  g(:, end+1, :) = 0;
