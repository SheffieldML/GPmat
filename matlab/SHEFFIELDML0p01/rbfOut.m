function [Y, G] = rbfOut(model, X);

% RBFOUT Output of an RBF model.
%
%	Description:
%
%	Y = RBFOUT(MODEL, X) gives the output of a radial basis function
%	model, the function is a wrapper for rbffwd.
%	 Returns:
%	  Y - the output.
%	 Arguments:
%	  MODEL - the model for which the output is required.
%	  X - the input data for which the output is required.
%
%	[Y, G] = RBFOUT(MODEL, X) gives the output of a radial basis
%	function model.
%	 Returns:
%	  Y - the output.
%	  G - the hidden layer activations.
%	 Arguments:
%	  MODEL - the model for which the output is required.
%	  X - the input data for which the output is required.
%	
%
%	See also
%	RBFFWD, RBF, MODELOUT


%	Copyright (c) 2006, 2007, 2008 Neil D. Lawrence


  if nargout > 1
    [Y, G] = rbffwd(model, X);
  else
    Y = rbffwd(model, X);
  end
