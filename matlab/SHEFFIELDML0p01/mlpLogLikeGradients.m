function g = mlpLogLikeGradients(model)

% MLPLOGLIKEGRADIENTS Multi-layer perceptron gradients.
%
%	Description:
%
%	G = MLPLOGLIKEGRADIENTS(MODEL) computes the gradients of the log
%	likelihood of a multi-layer perceptron with respect to the
%	parameters. For networks with one hidden layer this is done by
%	wrapping the mlpgrad command.
%	 Returns:
%	  G - the gradients of the model log likelihood.
%	 Arguments:
%	  MODEL - the model structure for computing the log likelihood.
%	
%
%	See also
%	MODELLOGLIKEIHOOD, MLPGRAD


%	Copyright (c) 2006 Neil D. Lawrence


if length(model.hiddenDim) == 1
  g = -mlpgrad(model, model.X, model.y);
else
  error('Not yet implemented.')
end