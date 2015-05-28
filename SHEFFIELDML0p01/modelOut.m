function [Y, Phi] = modelOut(model, X, varargin)

% MODELOUT Give the output of a model for given X.
%
%	Description:
%
%	Y = MODELOUT(MODEL, X) gives the output of the model for a given
%	input X. For latent variable models it gives a position in data
%	space given a position in latent space.
%	 Returns:
%	  Y - output location(s) corresponding to given input locations.
%	 Arguments:
%	  MODEL - structure specifying the model.
%	  X - input location(s) for which output is to be computed.
%
%	[PHI, Y] = MODELOUT(MODEL, X) gives the output of the model for a
%	given input X. For latent variable models it gives a position in
%	data space given a position in latent space.
%	 Returns:
%	  PHI - output basis function(s) corresponding to given input
%	  Y - output location(s) corresponding to given input locations.
%	 Arguments:
%	  MODEL - structure specifying the model.
%	  X - input location(s) for which output is to be computed.
%	
%	
%
%	See also
%	MODELCREATE


%	Copyright (c) 2005, 2006 Neil D. Lawrence


%	With modifications by Cark Henrik Ek 2008


fhandle = str2func([model.type 'Out']);
if nargout > 1
  [Y, Phi] = fhandle(model, X, varargin{:});
else
  Y = fhandle(model, X, varargin{:});
end
if(isfield(model,'indexOut')&&~isempty(model.indexOut))
  Y(:,setdiff(1:1:size(Y,2),model.indexOut)) = NaN;
end
