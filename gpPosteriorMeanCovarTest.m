function [analMu, analCov, diffMu, diffCov] = gpPosteriorMeanCovarTest(model, X)
% GPPOSTERIORMEANCOVARTEST Test the gradients of the mean and covariance.
% FORMAT
% DESC tests the gpPosteriorMeanCovar and gpPosteriorGradMeanCovar
% functions.
% ARG model : the model to test the gradients for.
% ARG X : the input locations to test the gradients for.
% RETURN analMu : the analytical gradients of the mean with respect
% to X.
% RETURN analCov : the analytical gradients of the covariance with respect
% to X.
% RETURN diffMu : the numerical gradients of the mean with respect
% to X.
% RETURN diffCov : the numerical gradients of the covariance with respect
% to X.
% RETURN delta : erros between the numerical and analytical
% gradients.
%
% SEEALSO : gpPosteriorMeanCovar, gpPosteriorGradMeanCovar
%
% COPYRIGHT : Neil D. Lawrence, 2006

% GP

[analMu, analCov] = gpPosteriorGradMeanCovar(model, X);
origX = X;
change = 1e-6;
for i = 1:size(X, 1)
  for j = 1:size(X, 2)
    X(i, j) = origX(i, j) - change;
    [muMinus, covMinus] = gpPosteriorMeanCovar(model, X);
    X(i, j) = origX(i, j) + change;
    [muPlus, covPlus] = gpPosteriorMeanCovar(model, X);
    X(i, j) = origX(i, j);
    diffMu{j}(i, :) = (muPlus(i, :) - muMinus(i, :))/(2*change);
    for k = 1:model.d
      % Not sure why the transpose is need here ... hope it isn't
      % two wrongs making a right ...
      diffCov{j, k}(:, i) = (covPlus{k}(i, :) - covMinus{k}(i, :))'/(2*change);
    end
    
  end
end
