% DROSPLOTCHIPDISTANCES Plot the ChIP binding site distance figures appearing in the paper
% FORMAT
% DESC Plot the ChIP binding site distance figures appearing in the paper
%
% COPYRIGHT : Antti Honkela, 2009

% SHEFFIELDML

drosLoadData;
demRunRankings;

%clear r;
t = [100, 250];
tfnames = {'twi', 'mef2'};

defaultorder = get(0, 'DefaultAxesColorOrder');
v = colormap;
%set(0, 'DefaultAxesColorOrder', [v([1, 22, 43, 64], :); 0, 0, 0]);
set(0, 'DefaultAxesColorOrder', [v([1, 17, 33, 48, 64], :); 0, 0, 0]);
set(0,'DefaultAxesLineStyleOrder',{'-','-','-','-','-',':'})

if ~exist('r'),
  for k=1:2,
    for l=1:2,
      tf = tfnames{l};
      r(k).(tf) = drosComputeChipDistanceCurve(chipdistances, drosexp, drosinsitu, {indrank.(tf), disimrank.(tf), oderank.(tf), mutarank.(tf), corrrank.(tf)}, t(k), tf);
    end
  end
end

if ~exist('r_focus'),
  for k=1:2,
    for l=1:2,
      tf = tfnames{l};
      r_focus(k).(tf) = drosComputeChipDistanceCurve(chipdistances, drosexp, drosinsitu, {indrank.(tf), disimrank.(tf), oderank.(tf), mutarank.(tf), corrrank.(tf)}, t(k), tf, 1);
    end
  end
end

FONTSIZE = 7;
rs = {r, r_focus};
labels = {'Global', 'Focused'};

figure(1);
j=1;
k=2;
for l=1:2,
  tf = tfnames{l};
  subplot(2, 2, l);
  semilogx(10.^[1:.1:6], 100*rs{j}(k).(tf))
  axis([10, 1e6, 0, 100])
  set(gca, 'FontSize', FONTSIZE);
  set(gca, 'XTick', 10.^[1:6]);
  title(sprintf('%s ChIP: %s, Top-%d', labels{j}, tf, t(k)))
  xlabel('Distance to CRM (bp)')
  if l==1,
    ylabel('Relative enrichment (%)')
  end
end
subplot(2, 2, 3:4);
plot(1:.1:6, 100*r(k).(tf))
axis off
set(gca, 'FontSize', FONTSIZE);
axis([10 11 0 1]);
legend('Single-target GP', 'Multiple-target GP', 'Single-target quadrature', 'Knock-outs', 'Correlation', 'Random', 'Location', 'North');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [20 20])
set(gcf, 'PaperPosition', [0 0 8.7 8.7])

figure(2);
for j=1:2,
  for k=1:2,
    for l=1:2,
      tf = tfnames{l};
      subplot(3, 4, 4*(k-1)+j+2*(l-1));
      semilogx(10.^[1:.1:6], 100*rs{j}(k).(tf))
      axis([10, 1e6, 0, 100])
      set(gca, 'FontSize', FONTSIZE);
      set(gca, 'XTick', 10.^[1:6]);
      title(sprintf('%s ChIP: %s, Top-%d', labels{j}, tf, t(k)))
      %xlabel('Genomic distance to closest ChIP binding site (bp)')
      xlabel('Distance to CRM (bp)')
      if l==1 && j==1,
	ylabel('Relative enrichment (%)')
      end
    end
  end
  subplot(3, 4, 10:11);
  plot(1:.1:6, 100*r(k).(tf))
  axis off
  set(gca, 'FontSize', FONTSIZE);
  axis([10 11 0 1]);
  legend('Single-target GP', 'Multiple-target GP', 'Single-target quadrature', 'Knock-outs', 'Correlation', 'Random', 'Location', 'North');
  set(gcf, 'PaperUnits', 'centimeters');
  set(gcf, 'PaperSize', [20 20])
  set(gcf, 'PaperPosition', [0 0 17.4 18])
end

set(0, 'DefaultAxesColorOrder', defaultorder);
set(0,'DefaultAxesLineStyleOrder',{'-'})
