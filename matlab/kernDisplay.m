function kernDisplay(kern)

% KERNDISPLAY Display the parameters of the kernel.

% IVM

feval([kern.type 'KernDisplay'], kern)