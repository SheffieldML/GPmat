function g = logLikeGradient(model, X, Y, varargin)

if(nargin<3)
  error('This function requires at least two arguments');
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

dlnZ_dmu = zeros(size(mu));

for j = 1:model.d
  dlnZ_dcovar{j} = zeros(size(X, 1));
  if missing
    ind = find(~isnan(Ydiff(:, j)));
    if length(ind) ~= 0
      U = jitChol(covarSigma(ind, ind));
    else
      U = [];
    end
  else
    ind = [1:size(Ydiff, 1)]';
  end
  UinvYdiff = U'\Ydiff(ind, j);
  dlnZ_dmu(ind, j) = U\UinvYdiff/factors(j);
  dlnZ_dcovar{j}(ind, ind) = 0.5*(UinvYdiff*UinvYdiff'/factors(j) - eye(length(ind)));
  dlnZ_dcovar{j}(ind, ind) = U\dlnZ_dcovar{j}(ind, ind)/U';
  dlnZ_dcovar{j}(ind, ind) = dlnZ_dcovar{j}(ind, ind)/factors(j);
end

g = zeros(size(X, 1), model.q);

for k = 1:model.d
  for j = 1:model.q
    g(:, j) = g(:, j) + dlnZ_dmu(:, k).*dmu{j}(:, k);
  end
  for j = 1:model.q
    % Since dlnZ_dcovar is symmetric we double dcovar apart from
    % the diagonal term ...
    diagCovar = diag(dcovar{j});
    dcovar2 = 2*dcovar{j};
    dcovar2(1:size(dcovar2, 2)+1:end) = diagCovar;
    for  i = 1:size(X, 1)
      g(i, j) = g(i, j) + dcovar2(:, i)'*dlnZ_dcovar{k}(:, i)*factors(k);
    end
  end
end

if isfield(model, 'dynamics') & ~isempty(model.dynamics)
    feval = str2func([model.dynamics.type,'SequenceLogLikeGradient']);
    if(isfield(model.dynamics,'indexIn')&&~isempty(model.dynamics.indexIn))
    dim = model.dynamics.indexIn;
  else
    dim = 1:1:size(g,2);
  end

  if isfield(model, 'balancing') & ~isempty(model.balancing)
    g(:,dim) = g(:,dim) + model.balancing*feval(model.dynamics, X);
  else
    g = g + feval(model.dynamics, X);
  end

	elseif isfield(model, 'prior') &  ~isempty(model.prior)
  
	for i = 1:size(X, 1)
		g(i, :) = g(i, :) + priorGradient(model.prior, X(i, :));
	end
end
   
g = (g(:)');
 