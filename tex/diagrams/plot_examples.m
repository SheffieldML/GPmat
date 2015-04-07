basetargets_twi2 = importdata('dros_twi_basetargets_public.txt');
m_gpsim_twi2 = drosGpsimLearn(drosexp, drosTF, 'twi', basetargets_twi2);
drosPlot2(m_gpsim_twi2, drosTF, drosexp)
set(gcf, 'PaperPosition', [0 0 8 6])
print -depsc2 gpsim_twi_example_basetargets_public
close all
m_gpdisim_twi2 = drosGpdisimLearn(drosexp, drosTF, 'twi', basetargets_twi2);
drosPlot2(m_gpdisim_twi2, drosTF, drosexp)
print -depsc2 gpdisim_twi_example_basetargets_public


samplegenes = {'FBgn0040600', 'FBgn0000927', 'FBgn0011206', 'FBgn0039286'};
r = load('dros_twi_gpsim_list_expro2_genes_dros_twi_basetargets_public_results.mat');
for k=1:4, II(k) = find(strcmp(r.targets, samplegenes{k})); end
l=3; m = drosGpsimRecreate(drosexp, drosTF, 'twi', [r.basetargets; r.targets{II(l)}], r.params{II(l)});
close all; drosPlot2(m, drosTF, drosexp)
set(gcf, 'PaperPosition', [0 0 6 4.5])
print -depsc2 ~/doc/masamb09/slides/diagrams/twi_gpsim_example3.eps
l=4; m = drosGpsimRecreate(drosexp, drosTF, 'twi', [r.basetargets; r.targets{II(l)}], r.params{II(l)});
close all; drosPlot2(m, drosTF, drosexp)
set(gcf, 'PaperPosition', [0 0 6 4.5])
print -depsc2 ~/doc/masamb09/slides/diagrams/twi_gpsim_example4.eps
