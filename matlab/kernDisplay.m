function kernDisplay(kern)

% KERNDISPLAY Display the parameters of the kernel.

% KERN

fhandle = str2func([kern.type 'KernDisplay']);
fhandle(kern);