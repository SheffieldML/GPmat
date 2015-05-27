function [res, accs] = drosBootstrapEvaluation(drosexp, drosinsitu, rankings, validation, repeats, threshold, do_isfilter)

% DROSBOOTSTRAPEVALUATION Performs bootstrap sampling of rankings
% FORMAT
% DESC Performs bootstrap sampling of rankings to assess significance of differences.
% ARG drosexp : the drosexp data structure from drosLoadData
% ARG drosinsitu : the drosinsitu data structure from drosLoadData
% ARG rankings : a cell array of R rankings, each of which is an array
% of indices of genes in drosexp in order of descending preference
% ARG validation : a binary vector of size(drosexp.genes) indicating
% the validation results of each gene, or NaN if the gene should be ignored
% ARG repeats : the number of repeats N
% ARG threshold : a vector of T n's in top-n to consider
% ARG is_filter : a flag whether to filter the results by positive in-situs
% (default: []=false)
% RETURN res : And RxRxT array where element (i,j,t) is the frequency that
% method i was better than method j for t'th threshold
% RETURN accs : An NxRxT array where element (n,i,t) is the accuracy of
% method i on repeat n for t'th threshold
%
% SEEALSO : drosLoadData, demRunRankings, drosDoBootstrap
%
% COPYRIGHT : Antti Honkela, 2009

% SHEFFIELDML

res = zeros([length(rankings), length(rankings), length(threshold)]);
accs = zeros([repeats, length(rankings), length(threshold)]);

baseids = find(~isnan(validation));

if do_isfilter,
  baseids = drosInsituPositives(baseids, drosinsitu, drosexp);
end

[B, I, J] = unique(drosexp.genes(baseids), 'first');

N = length(B);
val = validation(baseids(I));

for k=1:length(rankings),
  [myr, I2, J2] = intersect(baseids, rankings{k});
  [dummy, I3] = sort(J2);
  r{k} = J(I2(I3));
end

for j=1:length(threshold),
  for n=1:repeats,
    I = randint(N, 1, [1, N]);
    inds = accumarray(I, 1, [N, 1]);
    
    myacc = zeros(1, length(rankings));
    for k=1:length(rankings),
      J = cumsum(inds(r{k}));
      t = find(J >= threshold(j), 1);
      myinds = zeros(size(inds));
      myinds(r{k}(1:t-1)) = inds(r{k}(1:t-1));
      myinds(r{k}(t)) = inds(r{k}(t)) - J(t) + threshold(j);
      myacc(k) = sum(myinds .* val) ./ sum(myinds);
    end
    accs(n, :, j) = myacc;
  end

  for k=1:length(rankings),
    for l=k+1:length(rankings),
      res(k, l, j) = mean(accs(:, k, j) >= accs(:, l, j));
      res(l, k, j) = mean(accs(:, k, j) < accs(:, l, j));
    end
  end
end


function r = randint(rows, cols, range),

%r = ((range(1)):(range(2)))';
r = floor(range(2)*rand(rows, cols)) + range(1);
