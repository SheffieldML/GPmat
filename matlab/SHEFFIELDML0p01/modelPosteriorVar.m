function varsigma = modelPosteriorVar(model, X);

% MODELPOSTERIORVAR variances of the posterior at points given by X.
%
%	Description:
%
%	SIGMA = MODELPOSTERIORVAR(MODEL, X) returns the posterior  variance
%	for a given set of points.
%	 Returns:
%	  SIGMA - the variances of the posterior distributions.
%	 Arguments:
%	  MODEL - the model for which the posterior will be computed.
%	  X - the input positions for which the posterior will be computed.
%	
%
%	See also
%	MODELCREATE, MODELPOSTERIORMEANVAR


%	Copyright (c) 2009 Neil D. Lawrence

  varExist = false;
  func = [model.type 'PosteriorVar'];
  if exist(func)==2
    varExist = true;
  end
  if ~varExist
    func = [model.type 'PosteriorMeanVar'];
  end
  fhandle = str2func(func);

  if str2num(version('-release'))>13
    if varExist
      varsigma = fhandle(model, X);
    else
      [mu, varsigma] = fhandle(model, X);
    end
  else 
    if varExist
      varsigma = feval(fhandle, model, X);
    else
      [mu, varsigma] = feval(fhandle, model, X);
    end
  end
end
  
