% DEMGPTWOSAMPLEEB Run GP two sample code on EB.

% GP

% Load data
load EB
t{1} = [0 times2 0 times2]';
y{1} = [c1mtEBA; c3mtEBA];
t{2} = [0 times2 0 times2]';
y{2} = [c1ptEBA; c3ptEBA];


[llr, models] = gpTwoSample(t, y);


[void, order] = sort(llr);
for i = order(end:-1:1)'
  fprintf([num2str(llr(i), 4) '\t' genes_vect{1, i} '\n'])
end

save EBresults.mat llr models 

