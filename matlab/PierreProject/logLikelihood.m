function ll = logLikelihood(model, X, Y, varargin)

if(nargin<3)
  error('This function requires at least two arguments.');
end

logTwoPi = log(2*pi);
if model.isMissingData
  [mu, covarSigma] = gpPosteriorMeanCovar(model, X);
  [dmu, dcovar] = gpPosteriorGradMeanCovar(model, X);
else
  [mu, covarSigma, factors] = gpPosteriorMeanCovar(model, X);
  [dmu, dcovar, factors] = gpPosteriorGradMeanCovar(model, X);
end

missing = true;
if ~any(isnan(Y))
	U = jitChol(covarSigma);
	missing = false;
end

Ydiff = (Y-mu);
ll =0;
for i = 1:model.d
if missing
  ind = find(~isnan(Ydiff(:, i)));
  if length(ind) ~= 0    
	U = jitChol(covarSigma(ind, ind));
  end
else
  ind = [1:size(Ydiff, 1)]';
end
if length(ind) ~= 0    
  UinvYdiff = U'\Ydiff(ind, i);
  logDet = logdet([], U);
  ll = ll + logDet + log(factors(i)) + (UinvYdiff'*UinvYdiff)/factors(i);
  ll = ll + logTwoPi*length(ind);
end
end

ll = -0.5*ll;
    
if isfield(model, 'dynamics') & ~isempty(model.dynamics)
  
	feval = str2func([model.dynamics.type 'SequenceLogLikelihood']);
	if isfield(model, 'dynamicsBalancing') & ~isempty(model.dynamicsBalancing)
    ll = ll + model.dynamicsBalancing*feval(model.dynamics, X, varargin{:});
  else
    ll = ll + feval(model.dynamics, X, varargin{:});
  end

	elseif isfield(model, 'prior') &  ~isempty(model.prior)
  
	for i = 1:size(X, 1)
		ll = ll + priorLogProb(model.prior, X(i, :));
	end
end


