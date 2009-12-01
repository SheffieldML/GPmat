function drosPrintBootstrapMatrices(m, values, ignore),

if nargin < 3,
  ignore = [];
end

N = size(m, 3);
SIZES = [20, 100, 250];

for k=1:N,
  print_it(m(values, values, k), SIZES(k), ignore);
end


function print_it(m, mysize, ignore),

[M, N] = size(m);

LABELS = {'ST', 'MT', 'CO', 'KO'};

fprintf('\\subfloat[][Top %d]{\n', mysize);
fprintf('\\begin{tabular}{lcccc}\n');
fprintf('    & ST  & MT  & CO  & KO \\\\\n');

for k=1:M,
  fprintf('%s', LABELS{k});
  for l=1:N,
    if any(k == ignore) || any(l == ignore),
      fprintf(' & -  ');
    elseif m(k,l) > .99,
      fprintf(' & ***');
    elseif m(k,l) > .95,
      fprintf(' & ** ');
    elseif m(k,l) > .9,
      fprintf(' & *  ');
    elseif m(k,l) > .8,
      fprintf(' & +  ');
    elseif m(k,l) > .7,
      fprintf(' & .  ');
    else
      fprintf(' &    ');
    end
  end
  fprintf('\\\\ \n')
end

fprintf('\\end{tabular}}\n');




% drosPrintBootstrapMatrices(matrices{1}, 1:4)
% drosPrintBootstrapMatrices(matrices{2}, 1:4, 4)
% drosPrintBootstrapMatrices(matrices{3}, 1:4)
% drosPrintBootstrapMatrices(matrices{4}, 1:4)
% drosPrintBootstrapMatrices(matrices{5}, 1:4, 4)
% drosPrintBootstrapMatrices(matrices{6}, 1:4)
