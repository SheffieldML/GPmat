function [mu, varsigma] = modelPosteriorMeanVar(model, X);

% MODELPOSTERIORMEANVAR Mean and variances of the posterior at points given by X.
%
%	Description:
%
%	[MU, SIGMA] = MODELPOSTERIORMEANVAR(MODEL, X) returns the posterior
%	mean and variance for a given set of points.
%	 Returns:
%	  MU - the mean of the posterior distribution.
%	  SIGMA - the variances of the posterior distributions.
%	 Arguments:
%	  MODEL - the model for which the posterior will be computed.
%	  X - the input positions for which the posterior will be computed.
%	
%
%	See also
%	MODELCREATE


%	Copyright (c) 2009 Neil D. Lawrence


  fhandle = str2func([model.type 'PosteriorMeanVar']);
  if str2num(version('-release'))>13
    [mu, varsigma] = fhandle(model, X);
  else 
    [mu, varsigma] = feval(fhandle, model, X);
  end
  
  
