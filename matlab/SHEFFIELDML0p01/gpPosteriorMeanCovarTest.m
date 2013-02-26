function [analMu, analCov, diffMu, diffCov] = gpPosteriorMeanCovarTest(model, X)

% GPPOSTERIORMEANCOVARTEST Test the gradients of the mean and covariance.
%
%	Description:
%
%	[ANALMU, ANALCOV, DIFFMU, DIFFCOV, DELTA] =
%	GPPOSTERIORMEANCOVARTEST(MODEL, X) tests the gpPosteriorMeanCovar
%	and gpPosteriorGradMeanCovar functions.
%	 Returns:
%	  ANALMU - the analytical gradients of the mean with respect to X.
%	  ANALCOV - the analytical gradients of the covariance with respect
%	   to X.
%	  DIFFMU - the numerical gradients of the mean with respect to X.
%	  DIFFCOV - the numerical gradients of the covariance with respect
%	   to X.
%	  DELTA - erros between the numerical and analytical gradients.
%	 Arguments:
%	  MODEL - the model to test the gradients for.
%	  X - the input locations to test the gradients for.
%	
%
%	See also
%	GPPOSTERIORMEANCOVAR, GPPOSTERIORGRADMEANCOVAR


%	Copyright (c) 2006 Neil D. Lawrence


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