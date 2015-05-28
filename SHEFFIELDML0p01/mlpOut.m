function [Y, G, A] = mlpOut(model, X);

% MLPOUT Output of an MLP model.
%
%	Description:
%
%	Y = MLPOUT(MODEL, X) gives the output of a multi-layer perceptron
%	model, for single hidden layer models the function is a wrapper for
%	mlpfwd.
%	 Returns:
%	  Y - the output.
%	 Arguments:
%	  MODEL - the model for which the output is required.
%	  X - the input data for which the output is required.
%
%	[Y, G] = MLPOUT(MODEL, X) gives the output of a multi-layer
%	perceptron model.
%	 Returns:
%	  Y - the output.
%	  G - the hidden layer activations.
%	 Arguments:
%	  MODEL - the model for which the output is required.
%	  X - the input data for which the output is required.
%	
%
%	See also
%	MLPFWD, MLP, MODELOUT


%	Copyright (c) 2006, 2007 Neil D. Lawrence


if length(model.hiddenDim)==1
  if nargout > 1
    if nargout > 2
      [Y, G, A] = mlpfwd(model, X);
    else
      [Y, G] = mlpfwd(model, X);
    end
  else
    Y = mlpfwd(model, X);
  end
else
  ndata = size(x, 1);
  G{1} = tanh(X*model.w{1} + repmat(model.b{1}, numData, 1));
  A{1} = G{1}*model.w{2} + repmat(model.b{2}, numData, 1);
  for i = 2:length(model.numHidden)
    G{i} = tanh(A{i-1});
    A{i} = G{i}*model.w{i+1} + repmat(model.b{i+1}, numData, 1);
  end
  switch model.outfn
    
   case 'linear' 
    y = A{end};
    
   otherwise 
    error('Output function not implemented in multiple hidden layer model.')
    
  end
end    

  