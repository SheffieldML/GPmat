% DROSPLOTEVALUATION Plot the accuracy figures appearing in the paper
% FORMAT
% DESC Plot the accuracy figures appearing in the paper
%
% COPYRIGHT : Antti Honkela, 2009

% SHEFFIELDML

tfnames = {'twi', 'mef2'};
FONTSIZE = 7;
styles = {'bo-', 'rd--', 'm*--', 'gs--', 'k--'};
t = [20, 100, 250];

rankings = {};
for k=1:length(tfnames),
  tf = tfnames{k};
  if exist('INCLUDE_TSNI') && INCLUDE_TSNI,
    rankings{k} = {indrank.(tf), disimrank.(tf), tsnirank.(tf), mutarank.(tf), corrrank.(tf)};
    nonmutarankings{k} = {indrank.(tf), disimrank.(tf), tsnirank.(tf), [], corrrank.(tf)};
    INCLUDE_EXTRA=1;
    EXTRA_LABEL={'Multiple-target GP', 'TSNI'};
  elseif exist('INCLUDE_ODE') && INCLUDE_ODE,
    rankings{k} = {indrank.(tf), oderank.(tf), tsnirank.(tf), mutarank.(tf), corrrank.(tf)};
    nonmutarankings{k} = {indrank.(tf), oderank.(tf), tsnirank.(tf), [], corrrank.(tf)};
    INCLUDE_EXTRA=1;
    EXTRA_LABEL={'Single-target quadrature', 'TSNI'};
  else
    rankings{k} = {indrank.(tf), disimrank.(tf), oderank.(tf), mutarank.(tf), corrrank.(tf)};
    nonmutarankings{k} = {indrank.(tf), disimrank.(tf), oderank.(tf), [], corrrank.(tf)};
    INCLUDE_EXTRA=0;
    EXTRA_LABEL='';
  end
end

clear accs;
clear pvals;
figure(1);
for k=1:2,
  if INCLUDE_EXTRA,
    subplot(3, 4, 2*k-1);
  else
    subplot(3, 5, 2*k-1);
  end
  tf = tfnames{k};

  [accs(:, :, k), pvals(:, :, k)] = drosPlotAccuracyBars(...
      rankings{k}, ...
      chip_validation.(tf), t, styles, [], drosexp, drosinsitu);
  set(gca, 'FontSize', FONTSIZE);
  title(sprintf('Global ChIP: %s', tfnames{k}), 'FontSize', FONTSIZE);
  xlabel('Top N to consider', 'FontSize', FONTSIZE);
  if k==1,
    ylabel('Relative enrichment (%)', 'FontSize', FONTSIZE);
  end
end

clear accs;
clear pvals;
for k=1:2,
  if INCLUDE_EXTRA,
    subplot(3, 4, 3+2*k);
  else
    subplot(3, 5, 4+2*k);
  end
  tf = tfnames{k};

  [accs(:, :, k), pvals(:, :, k)] = drosPlotAccuracyBars(...
      nonmutarankings{k}, ...
      mutant_validation.(tf), t, styles, [], drosexp, drosinsitu);
  set(gca, 'FontSize', FONTSIZE);
  title(sprintf('Global knock-outs: %s', tfnames{k}), 'FontSize', FONTSIZE);
  xlabel('Top N to consider', 'FontSize', FONTSIZE);
  if k==1,
    ylabel('Relative enrichment (%)', 'FontSize', FONTSIZE);
  end
end

clear accs;
clear pvals;
for k=1:2,
  if INCLUDE_EXTRA,
    subplot(3, 4, 2*k);
  else
    subplot(3, 5, 2*k);
  end
  tf = tfnames{k};

  [accs(:, :, k), pvals(:, :, k)] = drosPlotAccuracyBars(...
      rankings{k}, ...
      chip_validation.(tf), t, styles, 1, drosexp, drosinsitu);
  set(gca, 'FontSize', FONTSIZE);
  title(sprintf('Focused ChIP: %s', tfnames{k}), 'FontSize', FONTSIZE);
  xlabel('Top N to consider', 'FontSize', FONTSIZE);
end

clear accs;
clear pvals;
for k=1:2,
  if INCLUDE_EXTRA,
    subplot(3, 4, 4+2*k);
  else
    subplot(3, 5, 5+2*k);
  end
  tf = tfnames{k};

  [accs(:, :, k), pvals(:, :, k)] = drosPlotAccuracyBars(...
      nonmutarankings{k}, ...
      mutant_validation.(tf), t, styles, 1, drosexp, drosinsitu);
  set(gca, 'FontSize', FONTSIZE);
  title(sprintf('Focused knock-outs: %s', tfnames{k}), 'FontSize', FONTSIZE);
  xlabel('Top N to consider', 'FontSize', FONTSIZE);
end


if INCLUDE_EXTRA,
  subplot(3, 4, [10, 11]);
  bar(rand(length(rankings{1})));
else
  subplot(3, 5, [5, 10]);
  bar(rand(length(rankings{1})));
end
hold on
plot(1:2, 1:2, 'k-.');
plot(1:2, 1:2, 'k--');

axis([-10 -9 -10 -9]);
set(gca, 'FontSize', FONTSIZE);
axis off
if INCLUDE_EXTRA,
  legend('Single-target GP', EXTRA_LABEL{:}, ...
	 'Knock-outs', 'Correlation', 'Filtered', 'Random', 'Location', 'North');
else
  legend(sprintf('Single-target GP'), sprintf('Multiple-target GP'), ...
	 sprintf('Single-target\nquadrature'), 'Knock-outs', ...
	 'Correlation', 'Filtered', 'Random', 'Location', 'West');
end
set(gca, 'FontSize', FONTSIZE);

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [20 20])
set(gcf, 'PaperPosition', [0 0 20 18])
hold off



if ~INCLUDE_EXTRA,
  figure(2);
  tf = 'twi';

  clear accs;
  clear pvals;

  rankings = {indrank.(tf), ssrank.(tf)};

  subplot(2, 3, 1);
  [accs(:, :, k), pvals(:, :, k)] = drosPlotAccuracyBars(rankings, ...
    chip_validation.(tf), t, styles, [], drosexp, drosinsitu);
  set(gca, 'FontSize', FONTSIZE);
  title(sprintf('Global ChIP: %s', tf));
  xlabel('Top N to consider');
  ylabel('Relative enrichment (%)');

  subplot(2, 3, 2);
  [accs(:, :, k), pvals(:, :, k)] = drosPlotAccuracyBars(rankings, ...
    mutant_validation.(tf), t, styles, [], drosexp, drosinsitu);
  set(gca, 'FontSize', FONTSIZE);
  title(sprintf('Global knock-outs: %s', tf));
  xlabel('Top N to consider');
  %ylabel('Relative enrichment (%)');

  subplot(2, 3, 4);
  [accs(:, :, k), pvals(:, :, k)] = drosPlotAccuracyBars(rankings, ...
    chip_validation.(tf), t, styles, 1, drosexp, drosinsitu);
  set(gca, 'FontSize', FONTSIZE);
  title(sprintf('Focused ChIP: %s', tf));
  xlabel('Top N to consider');
  ylabel('Relative enrichment (%)');

  subplot(2, 3, 5);
  [accs(:, :, k), pvals(:, :, k)] = drosPlotAccuracyBars(rankings, ...
    mutant_validation.(tf), t, styles, 1, drosexp, drosinsitu);
  set(gca, 'FontSize', FONTSIZE);
  title(sprintf('Focused knock-outs: %s', tf));
  xlabel('Top N to consider');
  %ylabel('Relative enrichment (%)');

  subplot(2, 3, [3,6]);
  bar(rand(length(rankings)));
  hold on
  plot(1:2, t(1:2), 'k--');
  plot(1:2, 1:2, 'k-.');

  axis([-10 -9 -10 -9]);
  set(gca, 'FontSize', FONTSIZE);
  axis off
  legend('n=12', 'n=7', ...
	 'Random', 'Filtered', 'Location', 'West');

  set(gcf, 'PaperUnits', 'centimeters');
  set(gcf, 'PaperSize', [20 20])
  set(gcf, 'PaperPosition', [0 0 9.65 9.65])
  %set(gcf, 'PaperSize', [20 20])
  %set(gcf, 'PaperPosition', [0 0 8.7 8.7])
  hold off
end
