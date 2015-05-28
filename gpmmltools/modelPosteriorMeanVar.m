function [mu, varsigma] = modelPosteriorMeanVar(model, X);

% MODELPOSTERIORMEANVAR Mean and variances of the posterior at points given by X.
% FORMAT
% DESC returns the posterior mean and variance for a given set of
% points.
% ARG model : the model for which the posterior will be computed.
% ARG x : the input positions for which the posterior will be
% computed.
% RETURN mu : the mean of the posterior distribution.
% RETURN sigma : the variances of the posterior distributions.
%
% SEEALSO : modelCreate
%
% COPYRIGHT : Neil D. Lawrence, 2009

% MLTOOLS

  fhandle = str2func([model.type 'PosteriorMeanVar']);
  if str2num(version('-release'))>13
    [mu, varsigma] = fhandle(model, X);
  else 
    [mu, varsigma] = feval(fhandle, model, X);
  end
  
  
