function kernDisplay(kern)

% KERNDISPLAY Display the parameters of the kernel.

% KERN


feval([kern.type 'KernDisplay'], kern)