function [Y, Phi] = ppcaOut(model, X);

% PPCAOUT Output of an PPCA model.
%
%	Description:
%
%	Y = PPCAOUT(MODEL, X) gives the output of a density network for a
%	given input.
%	 Returns:
%	  Y - the output.
%	 Arguments:
%	  MODEL - the model for which the output is required.
%	  X - the input data for which the output is required.
%	
%
%	See also
%	PPCACREATE, MODELOUT


%	Copyright (c) 2008 Neil D. Lawrence



if nargin<2  
  % If we are just updating outer layer, basis functions are stored.
  Y = model.X*model.W + repmat(model.b, model.M, 1);
  if nargout>1
    Phi = model.X;
  end
else
  Y = X*model.W + repmat(model.b, size(X, 1), 1);
  if nargout>1
    Phi = X;
  end
end