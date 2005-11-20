function k = kernDiagCompute(kern, x)

% KERNELCOMPUTE Compute the kernel given the parameters and X.

% KERN

fhandle = str2func([kern.type 'KernDiagCompute']);
k = fhandle(kern, x);
