function g = gpLogLikeGradX(model)

% GPLOGLIKEGRADX Gradient of the log likelihood with respect to X.

% GP

K = kernCompute(model.kern, model.X);
xDim = size(model.X, 2);
numData = size(model.X, 1);
g = zeros(numData, xDim);
fhandle = str2func([model.type 'CovarianceGradient']);
% gx can take up a lot of memory
gx = gpXGradient(model);

for i = 1:size(model.m, 2)
  if any(isnan(model.m(:, i)));
    ind = find(isnan(model.m(:, i)));
    tm = model.m(:, i);
    tm(ind, :) = [];
    tK = K;
    tK(ind, :) = [];
    tK(:, ind) = [];
    tinvK = pdinv(tK);
    covGrad = fhandle(tinvK, tm);
    counter = 0;
    for j = 1:numData
      if ~isnan(model.m(j, i))
        counter = counter+1;
        for k = 1:xDim
          tgx = gx(:, k, j);
          tgx(ind) = [];
          g(j, k) = g(j, k)+sum(tgx.*covGrad(:, counter));
        end
      end
    end
  else
    invK = pdinv(K);
    covGrad = fhandle(invK, model.m(:,i));
    for j = 1:numData
      for k = 1:xDim
        g(j, k) = g(j, k)+gx(:, k, j)'*covGrad(:, j);
      end
    end  
  end  
end
g = g(:)';