function [Y, Phi] = dnetOut(model, X);

% DNETOUT Output of an DNET model.
%
%	Description:
%
%	Y = DNETOUT(MODEL, X) gives the output of a density network for a
%	given input.
%	 Returns:
%	  Y - the output.
%	 Arguments:
%	  MODEL - the model for which the output is required.
%	  X - the input data for which the output is required.
%	
%
%	See also
%	DNETCREATE, MODELOUT


%	Copyright (c) 2008 Neil D. Lawrence


if nargin< 2  
  % Assume projection of latent samples is required (model.X_u).
  if model.basisStored && isfield(model, 'Phi')
    % If we are just updating outer layer, basis functions are stored.
    Y = model.Phi*model.A + repmat(model.b, model.M, 1);
    if nargout>1
      Phi = model.Phi;
    end
  else
    if nargout < 2
      Y = modelOut(model.mapping, model.X_u);
    else
      [Y, Phi] = modelOut(model.mapping, model.X_u);
    end
  end
else
  if nargout < 2
    Y = modelOut(model.mapping, X);
  else
    [Y, Phi] = modelOut(model.mapping, X);
  end
end