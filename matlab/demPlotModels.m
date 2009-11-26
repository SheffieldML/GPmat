% DEMPLOTMODELS Plot sample model figures appearing in the paper
% FORMAT
% DESC Plot sample model figures appearing in the paper
%
% COPYRIGHT : Antti Honkela, 2009

% DISIMRANK

drosLoadData;

tf = 'twi';

plottargets0 = {'FBgn0003486', 'FBgn0033188', 'FBgn0035257'};
plottargets = {};
close all;
for k=1:length(plottargets0),
  plottargets{k} = drosexp.probes{strcmp(drosexp.genes, plottargets0{k})};

  m = drosGpdisimLearn(drosexp, drosTF, tf, plottargets{k});
  drosPlot(m, length(plottargets0), k)
end
figure(3); 
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperPosition', [0, 0, 10, 12])
% print -depsc2 dros_gpdisim_singletarget_models
fprintf('Press any key to continue...');
pause;
if ~exist('m_multiplot'),
  m_multiplot = drosGpdisimLearn(drosexp, drosTF, tf, plottargets);
end
close all; drosPlot(m_multiplot)
figure(3); 
set(gcf, 'PaperUnits', 'centimeters')
set(gcf, 'PaperPosition', [0, 0, 10, 12])
% print -depsc2 dros_gpdisim_multitarget_model
