function kx = kerneldiag(x, theta)

% KERNELDIAG Compute the diagonal of the kernel function

kx = ones(size(x, 1), 1)*theta(2);
