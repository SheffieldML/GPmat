% DEMGPTWOSAMPLELIF Run GP two sample code on LIF.

% GP

% Load data
load LIF
t{1} = [times1 times1]';
y{1} = [c1mtLifA; c3mtLifA];
t{2} = [times1 times1]';
y{2} = [c1ptLifA; c3ptLifA];

[llr, models] = gpTwoSample(t, y);

[void, order] = sort(llr);
for i = order(end:-1:1)'
  fprintf([num2str(llr(i), 4) '\t' genes_vect{1, i} '\n'])
end

save LIFresults.mat llr models 
