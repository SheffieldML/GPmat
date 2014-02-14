% DEMPLOTMEF2MODELS Plot Mef2 sample model figures appearing in the paper
% FORMAT
% DESC Plot Mef2 sample model figures appearing in the paper
%
% COPYRIGHT : Antti Honkela, 2009

% SHEFFIELDML

drosLoadData;

tf = 'mef2';

plottargets0 = {'FBgn0010434', 'FBgn0030955', 'FBgn0035767'};
plottargets = {};
for k=1:length(plottargets0),
  plottargets{k} = drosexp.probes{strcmp(drosexp.genes, plottargets0{k})};
end
close all;
if ~exist('m'),
  for k=1:length(plottargets),
    m{k} = drosGpdisimLearn(drosexp, drosTF, tf, plottargets{k});
  end
end
for k=1:length(plottargets);
  drosPlot(m{k}, length(plottargets), k)
end
figure(3); 
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperPosition', [0, 0, 10, 6])
% print -depsc2 dros_gpdisim_singletarget_models
