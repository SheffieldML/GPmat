function gX = computeGradKernel(x, lntheta, type, x2)

% COMPUTEGRADKERNEL Compute the gradient of kernel given wrt X.


%/~ This function should return something with same number of columns as
%X and the same number of rows as x2 which means its the size of x2
%~/
if size(x, 1) > 1
  error('x should be a row vector')
end
lntheta=log(thetaConstrain(exp(lntheta)));
theta = exp(lntheta);
gX = zeros(size(x2));
switch type
 case 'linear'
  linearGrad = x2*theta(1);
  gX = linearGrad;  
  
 case 'rbf'
  rbfGrad = zeros(size(x2));
  n2 = dist2(x2, x);
  wi2 = (.5 .* theta(1));
  rbfPart = theta(2)*exp(-n2*wi2);
  for i = 1:size(x, 2)
    rbfGrad(:, i) = theta(1)*(x2(:, i) - x(i)).*rbfPart;
  end
  gX = rbfGrad;
  
 case 'regular'
  rbfGrad = zeros(size(x2));
  n2 = dist2(x2, x);
  wi2 = (.5 .* theta(1));
  rbfPart = theta(2)*exp(-n2*wi2);
  for i = 1:size(x, 2)
    rbfGrad(:, i) = theta(1)*(x2(:, i) - x(i)).*rbfPart;
  end
  linearGrad = x2*theta(5);
  gX = rbfGrad + linearGrad;
  
 case 'ARD'
  rbfGrad = zeros(size(x2));
  scales = diag(sqrt(theta(6:(5+size(x, 2)))));
  x = x*scales;
  x2 = x2*scales;
  n2 = dist2(x2, x);
  wi2 = (.5 .* theta(1));
  rbfPart = theta(2)*exp(-n2*wi2);
  for i = 1:size(x, 2)
    rbfGrad(:, i) = theta(1)*scales(i, i)*(x2(:, i) - x(i)).*rbfPart;
  end
  linearGrad = x2*scales*theta(5);
  gX = rbfGrad + linearGrad;
end




