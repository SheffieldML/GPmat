function cmpndKernDisplay(kern)

% CMPNDKERNDISPLAY Display the parameters of the compound kernel.

% IVM

for i = 1:length(kern.comp)
  kernDisplay(kern.comp{i});
end