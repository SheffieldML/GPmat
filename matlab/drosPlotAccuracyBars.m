function [r, pvals] = drosPlotAccuracyBars(rankings, validation, t, styles, filter, drosexp, drosinsitu),

% DROSPLOTACCURACYBARS Plot accuracies of alternative ranking methods.
% FORMAT
% DESC Plot a bar chart of accuracies of different ranking methods
% with stars indicating the p-value of the result against a random baseline
% ARG rankings : A cell array of rankings, each of which is an array
% of indices of genes in drosexp in order of descending preference
% ARG validation : A binary vector of size(drosexp.genes) indicating
% the validation results of each gene, or NaN if the gene should be ignored
% ARG t : a list of numbers of top-ranking genes to consider
% ARG styles : a cell array of plotting styles (optional, set empty for default)
% ARG filter : a flag whether to filter the results by positive in-situs
% (default: []=false)
% ARG drosexp : drosexp data set returned by drosLoadData
% (only needed when filter is nonempty)
% ARG drosinsitu : drosinsitu data set returned by drosLoadData
% (only needed when filter is nonempty)
% RETURN r : accuracies of each ranking
% RETURN pvals : p-values of ranking results against a random baseline
%
% SEEALSO : drosLoadData
%
% COPYRIGHT : Antti Honkela, 2009

% SHEFFIELDML

if nargin < 4 || isempty(styles),
  styles = {'bo-', 'gx-', 'r+-', 'cd--', 'm*--', 'ks--', 'k--'};
  %styles = {'bo-', 'gx-', 'r+-', 'cd-', 'ks-', 'm*-', 'm^-', 'mv--', 'k--'};
  %styles = {'b', 'g', 'r', 'c', 'k', 'm', 'm--'};
end

if nargin < 5,
  filter = [];
  drosexp = [];
  drosinsitu = [];
end

if ~iscell(rankings),
  rankings = {rankings};
end

r = zeros(length(t), length(rankings));
pvals = zeros(length(t), length(rankings));
L = 0;

values = [];

val2 = validation;
val2(isnan(val2)) = 0;
if ~isempty(filter),
  M = length(drosRemoveDuplicateGenes(drosexp, drosInsituPositives(find(~isnan(validation)), drosinsitu, drosexp)));
  K = length(drosRemoveDuplicateGenes(drosexp, drosInsituPositives(find(val2), drosinsitu, drosexp)));
  M0 = length(drosRemoveDuplicateGenes(drosexp, find(~isnan(validation))));
  K0 = length(drosRemoveDuplicateGenes(drosexp, find(val2)));
  %M0 = sum(~isnan(validation));
  %K0 = nansum(validation);
else
  if nargin > 5,
    M = length(drosRemoveDuplicateGenes(drosexp, find(~isnan(validation))));
    K = length(drosRemoveDuplicateGenes(drosexp, find(val2)));
  else
    M = sum(~isnan(validation));
    K = nansum(validation);
  end
end

for k=1:length(rankings),
  if ~isempty(filter),
    myranking = drosInsituPositives(rankings{k}, drosinsitu, drosexp);
  else
    myranking = rankings{k};
  end
  for l=1:length(t),
    if t(l) <= length(myranking),
      val = validation(myranking(1:t(l)));
      r(l, k) = nanmean(val);
      pvals(l, k) = 1 - hygecdf(nansum(val)-1, M, K, sum(~isnan(val)));
    else
      r(l, k) = 0;
      pvals(l, k) = 1;
    end
  end
end

%disp(values)
h = bar(100*r);
hold on
% Move the texts slightly because Matlab misaligns them
xfudge = -.05;
for k=1:length(rankings),
  x0 = get(get(h(k), 'Children'), 'XData');
  x = mean(x0(2:3, :)) + xfudge;
  y = get(h(k), 'YData');
  
  for l=1:length(x),
    if y(l) < 100 && y(l) > 0,
      if pvals(l, k) < .001,
	ht = text(x(l), y(l)+2, '***', 'Rotation', 90, 'VerticalAlignment', 'cap', 'HorizontalAlignment', 'left');
      elseif pvals(l, k) < .01,
	ht = text(x(l), y(l)+2, '**', 'Rotation', 90, 'VerticalAlignment', 'cap', 'HorizontalAlignment', 'left');
      elseif pvals(l, k) < .05,
	ht = text(x(l), y(l)+2, '*', 'Rotation', 90, 'VerticalAlignment', 'cap', 'HorizontalAlignment', 'left');
      end
    end
  end
end

%r(k+1) = b/d;

if ~isempty(filter),
  plot([0, length(t)+1], 100* K/M * [1 1], 'k-.');
  plot([0, length(t)+1], 100* K0/M0 * [1 1], 'k--');
else
  plot([0, length(t)+1], 100* K/M * [1 1], 'k--');
end
axis([.4, length(t)+.6, 0, 100])
hold off
set(gca, 'XTickLabel', t);
