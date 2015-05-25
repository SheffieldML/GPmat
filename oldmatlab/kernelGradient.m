function g = kernelGradient(parameter, x, h, A, type, prior)

if nargin < 6
  prior = 1;
end
[K, rbfPart, linearPart, dist2xx] = kernel(x, parameter, type);
invK = inv(K);
covGrad = covarianceGradient(invK, h, A);
g(1) = -.5*sum(sum(covGrad.*rbfPart.*dist2xx))*2*parameter(1);
g(2) = sum(sum(covGrad.*rbfPart/(parameter(2)*parameter(2))))*2*parameter(2);
g(3) = sum(sum(covGrad.*eye(size(x, 1))))*2*parameter(3);
g(4) = sum(sum(covGrad))*2*parameter(4);


switch type
 case 'regular'
  g(5) = sum(sum(covGrad.*(x*x')))*2*parameter(5);
 case 'ARD'
  scales = diag(parameter(6:(5+size(x, 2))));
  g(5) = sum(sum(covGrad.*(x*(scales*scales)*x')))*2*parameter(5);
  for i = 1:size(x, 2)
    g(5+i) = sum(sum(covGrad.*(parameter(5)*parameter(5)*x(:, i)*x(:, i)' ...
			       -.5*parameter(1)*parameter(1)*dist2(x(:, i), ...
				     x(:, i)) ...
			       .*rbfPart)))*2*scales(i, i);
  end
end
if prior
  g = g - 2./parameter;
end
g = -g;

