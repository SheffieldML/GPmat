function kernDisplay(kern)

% KERNDISPLAY Display the parameters of the kernel.

feval([kern.type 'KernDisplay'], kern)