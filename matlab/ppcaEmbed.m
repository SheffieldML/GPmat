function [X, sigma2, W] = ppcaEmbed(Y, dims)

% PPCAEMBED Embed data set with probabilistic PCA.
% FORMAT
% DESC returns latent positions for a given data set via probabilistic
% PCA.
% ARG Y : the data set which you want the latent positions for.
% ARG dims : the dimensionality of the latent space.
% RETURN X : the latent positions.
% RETURN sigma2 : the variance not explained by the latent positions.
% RETURN W : the matrix required to invert the transformation, Y=X*W'.
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% SEEALSO : lleEmbed, isomapEmbed

% MLTOOLS

if ~any(any(isnan(Y)))
  if size(Y, 1)<size(Y, 2)
    Ymean = mean(Y);
    Ycentre = zeros(size(Y));
    for i = 1:size(Y, 1)
      Ycentre(i, :) = Y(i, :) -Ymean;
    end
    if size(Ycentre, 2)>30000
      % Bug in MATLAB 7.0 means you have to do this.
      innerY = zeros(size(Ycentre, 1));
      for i = 1:size(Ycentre, 1)
        innerY(i, :) = Ycentre(i, :)*Ycentre';
      end
    else
      innerY = Ycentre*Ycentre';
    end
    [v, u] = eigdec(innerY, dims); 
    v(find(v<0))=0;
    X = u(:, 1:dims)*sqrt(size(Y, 1));
    v = v/sqrt(size(Y, 1));
    sigma2 = (trace(innerY) - sum(v))/(size(Y, 2)-dims);
    W = X'*Ycentre;

  else
    [v, u] = pca(Y);
    v(find(v<0))=0;
    Ymean = mean(Y);
    Ycentre = zeros(size(Y));
    for i = 1:size(Y, 2);
      Ycentre(:, i) = Y(:, i) - Ymean(i);
    end
    X = Ycentre*u(:, 1:dims)*diag(1./sqrt(v(1:dims)));
    sigma2 = mean(v(dims+1:end));
    W = X'*Ycentre;
  end
else
  % Hacky implementation of Probabilistic PCA for when there is missing data.
  iters = 100;
  % Initialise W randomly
  d = size(Y, 2);
  q = dims;
  N = size(Y, 1);
  W = randn(d, q)*1e-3;
  sigma2 = 1;
  mu = zeros(d, 1);
  for i = 1:d
    obs = find(~isnan(Y(:, i)));
    if length(obs)>0
      mu(i) = mean(Y(obs, i));
    else
      mu(i) = 0;
    end
  end
  numObs = sum(sum(~isnan(Y)));
  for i = 1:iters
    M = W'*W + sigma2*eye(q);
    invM = inv(M);
    exp_xxT = zeros(q);
    exp_x = zeros(N, q);
    for n = 1:N
      obs = find(~isnan(Y(n, :)));
      exp_x(n, :) = (invM*W(obs, :)'*(Y(n, obs)' - mu(obs)))';
    end
    exp_xxT = N*sigma2*invM + exp_x'*exp_x;
    s = zeros(d, q);
    s2 = 0;
    for n = 1:N
      obs = find(~isnan(Y(n, :)));
      subY = zeros(size(Y(n, :)))';
      subY(obs) = Y(n, obs)' - mu(obs);
      s = s + (subY)*exp_x(n, :);
      s2 = s2 + sum((Y(n, obs)' - mu(obs)).^2) - 2*exp_x(n, :)*W(obs, :)'*(Y(n, obs)' - mu(obs));
    end
    W = s*inv(exp_xxT);
    sigma2 = 1/(numObs)*(s2 + trace(exp_xxT*W'*W));
  end
  
  X = exp_x;
end
