function L = gplvmlikelihood(params, Y, X)

% GPLVMLIKELIHOOD Calculate the likelihood of the data set.

global SQRT % Tells us to optimise in sqrt of theta space

dataDim = size(Y, 2);
numData = size(Y, 1);

lntheta = params(end-2:end);
if nargin < 3
  latentDim = (length(params) -3)/numData;
  % We must have to optimise X as well -- it isn't given!
  X = reshape(params(1:numData*latentDim), numData, latentDim);
  xOpt = 1;
else
  xOpt = 0;
end
    
if SQRT
  theta = lntheta*lntheta;
else
  theta = exp(lntheta);
end

theta = thetaConstrain(theta);

[K, invK] = computeKernel(X, theta);
  
L = -dataDim*numData/2*log(2*pi)-dataDim/2*logdet(K);
for d= 1:dataDim
  L = L -0.5*Y(:, d)'*invK*Y(:, d);
end

% Include priors over parameters for MAP estimate

L = L - sum(log(theta));

if xOpt
  L = L - .5*sum(sum(X.*X));
end
L = -L;
