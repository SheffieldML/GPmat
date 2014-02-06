function drosPrintBootstrapMatrices(m, values, ignore),

% DROSPRINTBOOTSTRAPMATRICES Print bootstrap sampling results as a LaTeX table
% FORMAT
% DESC Print bootstrap sampling results as a LaTeX table
% ARG m : the res array returned by drosBoostrapEvaluation
% ARG values : the indices of rankings in m to consider
% ARG ignore : the indices of rankings in m to ignore, as indices within values
%
% SEEALSO : drosLoadData, demRunRankings, drosDoBootstrap, drosBootstrapEvaluation
%
% COPYRIGHT : Antti Honkela, 2009

% SHEFFIELDML

if nargin < 3,
  ignore = [];
end

N = size(m, 3);
SIZES = [20, 100, 250];

print_it(m(values, values, :), SIZES, ignore);


function print_it(m, mysize, ignore),

[M, N, O] = size(m);

LABELS = {'ST', 'MT', 'QR', 'CO', 'KO'};

fprintf('\\multicolumn{6}{c|}{Top %d} &\\multicolumn{6}{c|}{Top %d}  & \\multicolumn{6}{c}{Top %d}\\\\\n', mysize);
fprintf('    & ST  & MT  & QR  & CO  & KO &     & ST  & MT  & QR  & CO  & KO  &     & ST  & MT  & QR  & CO  & KO \\\\\n');

for k=1:M,
  for j=1:O,
    fprintf('%s', LABELS{k});
    for l=1:N,
      if any(k == ignore) || any(l == ignore),
	fprintf(' & -  ');
      elseif m(k,l,j) > .99,
	fprintf(' & ***');
      elseif m(k,l,j) > .95,
	fprintf(' & ** ');
      elseif m(k,l,j) > .9,
	fprintf(' & *  ');
      elseif m(k,l,j) > .8,
	fprintf(' & +  ');
      elseif m(k,l,j) > .7,
	fprintf(' & .  ');
      else
	fprintf(' &    ');
      end
    end
    if j<O,
      fprintf(' & ');
    end
  end
  fprintf('\\\\\n')
end




% drosPrintBootstrapMatrices(matrices{1}, 1:4)
% drosPrintBootstrapMatrices(matrices{2}, 1:4, 4)
% drosPrintBootstrapMatrices(matrices{3}, 1:4)
% drosPrintBootstrapMatrices(matrices{4}, 1:4)
% drosPrintBootstrapMatrices(matrices{5}, 1:4, 4)
% drosPrintBootstrapMatrices(matrices{6}, 1:4)
