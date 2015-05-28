% DROSPLOTCHIPDISTANCES Plot the accuracy figures appearing in the paper
% FORMAT
% DESC Plot the accuracy figures appearing in the paper
%
% COPYRIGHT : Antti Honkela, 2009

% SHEFFIELDML

tfnames = drosTF.names;
FONTSIZE = 8;
styles = {'bo-', 'rd--', 'm*--', 'gs--', 'k--'};
t = [20, 100, 250];

rankings = {};
for k=1:length(tfnames),
  tf = tfnames{k};
  nonmutarankings{k} = {indrank.(tf), [], corrrank.(tf)};
  if isfield(mutarank, tf),
    rankings{k} = {indrank.(tf), mutarank.(tf), corrrank.(tf)};
  else
    rankings{k} = nonmutarankings{k};
  end
end

clear accs;
clear pvals;
figure(1);
for k=1:length(tfnames),
  subplot(3, length(tfnames), k); %3*(k-1)+2);
  tf = tfnames{k};

  [accs(:, :, k), pvals(:, :, k)] = drosPlotAccuracyBars(...
      rankings{k}, ...
      chip_validation.(tf), t, styles, [], drosexp, drosinsitu);
  set(gca, 'FontSize', FONTSIZE);
  title(sprintf('Global ChIP: %s', tfnames{k}));
  xlabel('Top N to consider');
  if k==1,
    ylabel('Relative enrichment (%)');
  end
end

clear accs;
clear pvals;
for k=1:length(tfnames),
  subplot(3, length(tfnames), length(tfnames)+k);
  tf = tfnames{k};

  if isfield(mutant_validation, tf),
    [accs(:, :, k), pvals(:, :, k)] = drosPlotAccuracyBars(...
	nonmutarankings{k}, ...
	mutant_validation.(tf), t, styles, [], drosexp, drosinsitu);
    set(gca, 'FontSize', FONTSIZE);
    title(sprintf('Global knock-outs: %s', tfnames{k}));
    xlabel('Top N to consider');
    if k==1,
      ylabel('Relative enrichment (%)');
    end
  end
end


subplot(3, length(tfnames), 3*length(tfnames)-[1, 0]);
if exist('INCLUDE_TSNI') && INCLUDE_TSNI,
  bar(rand(5));
else
  bar(rand(length(rankings{1})));
end
hold on
plot(1:2, t(1:2), 'k--');

axis([-10 -9 -10 -9]);
set(gca, 'FontSize', FONTSIZE);
axis off
if exist('INCLUDE_TSNI') && INCLUDE_TSNI,
  legend('Single-target models', 'Multiple-target models', 'TSNI', ...
	 'Knock-outs', 'Correlation', 'Random', 'Location', 'North');
else
  legend('Single-target models', ... %'Multiple-target models', ...
	 'Knock-outs', 'Correlation', 'Random', 'Location', 'North');
end

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [20 20])
set(gcf, 'PaperPosition', [0 0 17.4 18])
hold off


clear accs;
clear pvals;
figure(2);
for k=1:length(tfnames),
  subplot(3, length(tfnames), k); %3*(k-1)+2);
  tf = tfnames{k};

  [accs(:, :, k), pvals(:, :, k)] = drosPlotAccuracyBars(...
      rankings{k}, ...
      chip_validation.(tf), t, styles, 1, drosexp, drosinsitu);
  set(gca, 'FontSize', FONTSIZE);
  title(sprintf('Focused ChIP: %s', tfnames{k}));
  xlabel('Top N to consider');
  if k==1,
    ylabel('Relative enrichment (%)');
  end
end

clear accs;
clear pvals;
for k=1:length(tfnames),
  subplot(3, length(tfnames), length(tfnames)+k);
  tf = tfnames{k};

  if isfield(mutant_validation, tf),
    [accs(:, :, k), pvals(:, :, k)] = drosPlotAccuracyBars(...
	nonmutarankings{k}, ...
	mutant_validation.(tf), t, styles, 1, drosexp, drosinsitu);
    set(gca, 'FontSize', FONTSIZE);
    title(sprintf('Focused knock-outs: %s', tfnames{k}));
    xlabel('Top N to consider');
    if k==1,
      ylabel('Relative enrichment (%)');
    end
  end
end


subplot(3, length(tfnames), 3*length(tfnames)-[1, 0]);
if exist('INCLUDE_TSNI') && INCLUDE_TSNI,
  bar(rand(5));
else
  bar(rand(length(rankings{1})));
end
hold on
plot(1:2, 1:2, 'k-.');
plot(1:2, 1:2, 'k--');

axis([-10 -9 -10 -9]);
set(gca, 'FontSize', FONTSIZE);
axis off
if exist('INCLUDE_TSNI') && INCLUDE_TSNI,
  legend('Single-target models', 'Multiple-target models', 'TSNI', ...
	 'Knock-outs', 'Correlation', 'Random', 'Location', 'North');
else
  legend('Single-target models', ... %'Multiple-target models', ...
	 'Knock-outs', 'Correlation', 'Random', 'Location', 'North');
end

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [20 20])
set(gcf, 'PaperPosition', [0 0 17.4 18])
hold off



% figure(3);
% tf = 'twi';

% clear accs;
% clear pvals;

% rankings = {indrank.(tf), ssrank.(tf)};

% subplot(2, 3, 1);
% [accs(:, :, k), pvals(:, :, k)] = drosPlotAccuracyBars(rankings, ...
%     chip_validation.(tf), t, styles, [], drosexp, drosinsitu);
% set(gca, 'FontSize', FONTSIZE);
% title(sprintf('Global ChIP: %s', tf));
% xlabel('Top N to consider');
% ylabel('Relative enrichment (%)');

% subplot(2, 3, 2);
% [accs(:, :, k), pvals(:, :, k)] = drosPlotAccuracyBars(rankings, ...
%     mutant_validation.(tf), t, styles, [], drosexp, drosinsitu);
% set(gca, 'FontSize', FONTSIZE);
% title(sprintf('Global knock-outs: %s', tf));
% xlabel('Top N to consider');
% %ylabel('Relative enrichment (%)');

% subplot(2, 3, 4);
% [accs(:, :, k), pvals(:, :, k)] = drosPlotAccuracyBars(rankings, ...
%     chip_validation.(tf), t, styles, 1, drosexp, drosinsitu);
% set(gca, 'FontSize', FONTSIZE);
% title(sprintf('Focused ChIP: %s', tf));
% xlabel('Top N to consider');
% ylabel('Relative enrichment (%)');

% subplot(2, 3, 5);
% [accs(:, :, k), pvals(:, :, k)] = drosPlotAccuracyBars(rankings, ...
%     mutant_validation.(tf), t, styles, 1, drosexp, drosinsitu);
% set(gca, 'FontSize', FONTSIZE);
% title(sprintf('Focused knock-outs: %s', tf));
% xlabel('Top N to consider');
% %ylabel('Relative enrichment (%)');

% subplot(2, 3, [3,6]);
% bar(rand(2));
% hold on
% plot(1:2, t(1:2), 'k--');
% plot(1:2, 1:2, 'k-.');

% axis([-10 -9 -10 -9]);
% set(gca, 'FontSize', FONTSIZE);
% axis off
% legend('n=12', 'n=7', ...
%        'Random', 'Filtered', 'Location', 'West');

% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperSize', [20 20])
% set(gcf, 'PaperPosition', [0 0 9.65 9.65])
% %set(gcf, 'PaperSize', [20 20])
% %set(gcf, 'PaperPosition', [0 0 8.7 8.7])
% hold off
