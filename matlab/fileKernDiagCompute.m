function k = fileKernDiagCompute(kern, x)

% FILEKERNDIAGCOMPUTE Load the diagonal of a kernel stored in a file.

% KERN

k = kern.variance*fileKernRead(kern, x, 'diag');
