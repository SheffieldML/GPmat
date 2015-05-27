function drosPlotExpPcts(drosexp, genes),
% DROSPLOTEXPPCTS Plot expression values
%
% Usage:
%   drosPlotExpPcts(drosexp, genes)
% COPYRIGHT : Antti Honkela, 2007
  
% SHEFFIELDML

if iscell(genes),
  genes = drosGetGeneinds(drosexp, genes);
end

N = length(genes);

pcts0 = drosexp.pctiles(genes, :, [2, 3, 4]);
labels = drosexp.genes(genes);

pcts = zeros(N, 12, 3);
pcts(:, :, 1) = min(reshape(pcts0(:, :, 1), [N, 12, 3]), [], 3);
pcts(:, :, 2) = mean(reshape(pcts0(:, :, 2), [N, 12, 3]), 3);
pcts(:, :, 3) = max(reshape(pcts0(:, :, 3), [N, 12, 3]), [], 3);

for k=1:N,
  if N>1,
    subplot(N, 1, k);
  end
  plot(1:12, exp(pcts(k, :, 2)));
  hold on;
  plot(1:12, exp(pcts(k, :, 1)), '--');
  plot(1:12, exp(pcts(k, :, 3)), '--');
  if k==N,
    set(gca, 'XTick', 1:12);
    set(gca, 'XTickLabel', drosexp.labels);
  else
    set(gca, 'XTick', 1:12);
    set(gca, 'XTickLabel', {});
  end
  axis tight
  ax = axis;
  axis([1, 12, ax(3:4)]);
  hold off
  text(0.5, mean(ax(3:4)), labels{k}, 'HorizontalAlignment', 'right');
end
