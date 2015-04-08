% IVMPUMADYNRESULTS

load demPumadyn1
dVals = [50 100 200 500 1000];
meanEr = mean(log10(testError));
stdEr = sqrt(var(log10(testError)));
a = errorbar(dVals, meanEr, stdEr);
set(a, 'linewidth', 2)
vals = 1:9;
ticks = log10([vals*1e-3 vals*1e-2 vals*1e-1 vals*1])

xLim = get(gca, 'xlim');
yLim = get(gca, 'ylim');

for i = 1:length(ticks)
  line(xLim, [ticks(i) ticks(i)], 'linewidth', 0.5, 'linestyle', ':')
end
xticks = get(gca, 'xtick')

for i = 1:length(xticks)
  line([xticks(i) xticks(i)], yLim,'linewidth', 0.5, 'linestyle', ':')
end
set(gca, 'ytick', ticks)
set(gca, 'xlim', xLim, 'ylim', yLim)
for i = 1:length(ticks)
  switch rem(i, 9)
   case {1, 2, 5}
     tickLabels{i} = num2str(10^ticks(i)) ;;
   otherwise
    tickLabels{i} = '';
  end
end
set(gca, 'yticklabel', tickLabels)
set(gca, 'fontname', 'times')
set(gca, 'fontsize', 18)

xlabel('d')
ylabel('sq. error')