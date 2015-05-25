function plot_exp(profs, legends),

STAGES = {'s2', 's4', 's5', 's8', 's9-10', ...
	  's11e', 's11l-12e', 's12e', 's12l', ...
	  's13e', 's13l', 's14'};

if isstruct(profs),
  N = size(profs.m, 1);

  v = 1 ./ mean(reshape(1 ./ profs.se.^2, [N, 12, 3]), 3);
  m = v .* mean(reshape(profs.m ./ (profs.se.^2), [N, 12, 3]), 3);
  v = sqrt(v);
else
  N = size(profs, 1);

  m = mean(reshape(profs, [size(profs, 1), 12, 3]), 3);
  v = std(reshape(profs, [size(profs, 1), 12, 3]), [], 3);
end

for k=1:N,
  subplot(N, 1, k);
  plot(1:12, m(k, :));
  hold on;
  plot(1:12, m(k, :)+v(k, :), '--');
  plot(1:12, m(k, :)-v(k, :), '--');
  if k==N,
    set(gca, 'XTick', 1:12);
    set(gca, 'XTickLabel', STAGES);
  else
    set(gca, 'XTick', 1:12);
    set(gca, 'XTickLabel', {});
  end
  axis tight
  ax = axis;
  axis([1, 12, ax(3:4)]);
  hold off
  if nargin > 1,
    text(0.5, mean(ax(3:4)), legends{k}, 'HorizontalAlignment', 'right');
  end
end
