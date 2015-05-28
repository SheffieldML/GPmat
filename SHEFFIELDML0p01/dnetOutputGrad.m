function g = dnetOutputGrad(model, X)

% DNETOUTPUTGRAD Evaluate derivatives of dnet model outputs with respect to parameters.
%
%	Description:
%
%	G = DNETOUTPUTGRAD(MODEL, X) evaluates the derivates of a density
%	network's outputs with respect to the parameters of the mapping
%	function.
%	 Returns:
%	  G - the gradient of the outputs of the density network perceptron
%	   with respect to each of the parameters. The size of the matrix is
%	   number of data x number of parameters x number of outputs of the
%	   model.
%	 Arguments:
%	  MODEL - the model for which the derivatives are to be computed.
%	  X - the input data locations where the gradients are to be
%	   computed.
%	
%
%	See also
%	MODELOUTPUTGRAD, DNETCREATE


%	Copyright (c) 2008 Neil D. Lawrence


  g = modelOutputGrad(model.mapping, X);
  % Add output gradient wrt beta.
  g(:, end+1, :) = 0;