% DATAPATH = sprintf('%s/data/eileen_data/magnus-data/', getenv('HOME'));
DATAPATH = sprintf('./data/eileen_data/');

expdata = importdata([DATAPATH 'expro/bdgp_expro_fbgn_mmgmos_refseq.txt']);
% expdata = importdata([DATAPATH, 'expro/bdgp_expro_fbgn_mmgmos_refseq.txt']);
expdata.data = expdata.data - repmat(mean(expdata.data), [size(expdata.data, 1), 1]) + mean(mean(expdata.data));

expse = importdata([DATAPATH, 'expro/bdgp_expro_fbgn_mmgmos_refseq_se.txt']);

prctiles = [5, 25, 50, 75, 95];
pcts = cell(5, 1);

for k=1:5,
  foo = importdata(sprintf('%s/expro/bdgp_expro_fbgn_mmgmos_refseq_prctile%d.txt', DATAPATH, prctiles(k)));
  pcts{k} = foo.data;
end

drosexp.pctvalues = prctiles;
drosexp.pctiles = cat(3, pcts{:});
drosexp.genes = foo.rowheaders;

drosexp.mean = expdata.data;
drosexp.se = expse.data;
drosexp.labels = {'s2', 's4', 's5', 's8', 's9-10', ...
          's11e', 's11l-12e', 's12e', 's12l', ...
          's13e', 's13l', 's14'};

chipdata = importdata([DATAPATH, 'chipchip/gene_binding_matrix.txt']);
droschip.labels = {'tin_2-4', 'tin_4-6', 'tin_6-8', 'bin_6-8', 'bin_8-10', 'bin_10-12', 'twi_2-4', 'twi_4-6', 'twi_6-8', 'bap_6-8', 'mef2_2-4', 'mef2_4-6', 'mef2_6-8', 'mef2_8-10', 'mef2_10-12'};
droschip.data = chipdata.data;
droschip.genes = chipdata.rowheaders;

drosTF.names = {'tin', 'bin', 'twi', 'bap', 'mef2'};
drosTF.labels = {'FBgn0004110', 'FBgn0045759', 'FBgn0003900', ...
		 'FBgn0004862', 'FBgn0011656'};
drosTF.chipinds = {1:3, 4:6, 7:9, 10, 11:15};

clear expdata expse pcts prctiles foo chipdata
