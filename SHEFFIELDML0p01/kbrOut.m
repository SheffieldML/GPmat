function [Y, G] = kbrOut(model, X);

% KBROUT Compute the output of a KBR model given the structure and input X.
%
%	Description:
%
%	Y = KBROUT(MODEL, X) computes the model parameters for the kernel
%	based regression model given inputs associated with rows and
%	columns.
%	 Returns:
%	  Y - the output results.
%	 Arguments:
%	  MODEL - the model structure for which the output is computed.
%	  X - the input data.
%
%	[Y, G] = KBROUT(MODEL, X) gives the output of a radial basis
%	function model.
%	 Returns:
%	  Y - the output.
%	  G - the values computed at the kernel.
%	 Arguments:
%	  MODEL - the model for which the output is required.
%	  X - the input data for which the output is required.
%	
%
%	See also
%	KBRCREATE, MODELCOMPUTE, MODELCREATE, KBREXPANDPARAM, KBREXTRACTPARAM


%	Copyright (c) 2005, 2006, 2008 Neil D. Lawrence


numData = size(X, 1);
if ~isfield(model, 'bias') & isfield(model, 'b')
  model.bias = model.b;
end
G = kernCompute(model.kern, X, model.X);
Y = G*model.A+ones(numData, 1)*model.bias;
