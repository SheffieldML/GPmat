% DEMGPTWOSAMPLEEBSH Run GP two sample code on EB data

% GP

% Load data
load EBsh

t{1} = [times_shEB times_shEB]';
y{1} = [c1shEB; c3shEB];
t{2} = [times_shEB times_shEB]';
y{2} = [a7shE13EB; c4shE13EB];


[llr, models] = gpTwoSample(t, y);


[void, order] = sort(llr);
for i = order(end:-1:1)'
  fprintf([num2str(llr(i), 4) '\t' genes_vect_shEB{1, i} '\n'])
end

save EBshresults.mat llr models 

