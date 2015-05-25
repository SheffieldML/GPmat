function kx = kernel(x, xKern, theta)

% KERNEL Compute the rbf kernel

n2 = dist2(x, xKern);
wi2 = theta(1)/2;
kx = theta(2)*exp(-n2*wi2);
   