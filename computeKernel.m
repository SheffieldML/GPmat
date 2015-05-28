function [kx, rbfPart, linearPart, n2] = computeKernel(x, lntheta, type, x2)

% COMPUTEKERNEL Compute the kernel given the parameters and X.

% KERN

% KERN

lntheta=log(thetaConstrain(exp(lntheta)));
theta = exp(lntheta);

if nargin < 3
  type = 'regular';
end

if nargin < 4
  
  switch type
   case 'linear'
    rbfPart = [];
    n2 = [];
    linearPart = x*x'*theta(1);
    kx = linearPart + theta(3)*eye(size(x, 1)) + theta(2);
    
   case 'rbf'
    n2 = dist2(x, x);
    wi2 = (.5 .* theta(1));
    rbfPart = theta(2)*exp(-n2*wi2);
    kx = rbfPart + theta(3)*eye(size(x, 1)) + theta(4);
    linearPart = [];
    
   case 'regular'
    n2 = dist2(x, x);
    wi2 = (.5 .* theta(1));
    rbfPart = theta(2)*exp(-n2*wi2);
    linearPart = x*x'*theta(5);
    kx = rbfPart + theta(3)*eye(size(x, 1)) + theta(4) + linearPart;
    
   case 'ARD'
    scales = diag(sqrt(theta(6:(5+size(x, 2)))));
    x = x*scales;
    n2 = dist2(x, x);
    wi2 = (.5 .* theta(1));
    rbfPart = theta(2)*exp(-n2*wi2);
    linearPart = x*x'*theta(5);
    kx = rbfPart + theta(3)*eye(size(x, 1)) + theta(4) + linearPart;
  end
  
else
  switch type
   case 'linear'
    rbfPart = [];
    n2 = [];
    linearPart = x*x2'*theta(1);
    kx = theta(2) + linearPart;
    
    
   case 'rbf'
    n2 = dist2(x, x2);
    wi2 = (.5 .* theta(1));
    rbfPart = theta(2)*exp(-n2*wi2);
    linearPart = [];
    kx = rbfPart + theta(4);
   
   case 'regular'
    n2 = dist2(x, x2);
    wi2 = (.5 .* theta(1));
    rbfPart = theta(2)*exp(-n2*wi2);
    linearPart = x*x2'*theta(5);
    kx = rbfPart + theta(4) + linearPart;
    
   case 'ARD'
    scales = diag(sqrt(theta(6:(5+size(x, 2)))));
    x = x*scales;
    x2 = x2*scales;
    n2 = dist2(x, x2);
    wi2 = (.5 .* theta(1));
    rbfPart = theta(2)*exp(-n2*wi2);
    linearPart = x*x2'*theta(5);
    kx = rbfPart + theta(4) + linearPart;
  end
end




