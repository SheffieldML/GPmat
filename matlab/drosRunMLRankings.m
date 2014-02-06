% DROSRUNMLRANKINGS Runs the baseline maximum likelihood ranking presented in the paper
% FORMAT
% DESC Runs the baseline maximum likelihood rankings presented in the
% paper for comparison
%
% COPYRIGHT : Antti Honkela, 2010

% SHEFFIELDML

drosLoadData
z_scores = drosexp.fitmean ./ sqrt(drosexp.fitvar);
active_genes = (mean(z_scores, 2) >= 1.8);

[r_twi, lls_twi, params_twi] = drosMLRank(drosexp, drosTF, 'twi', drosexp.probes(active_genes));
r = r_twi; lls = lls_twi; params = params_twi;
save results/dros_mlode_twi_list_expro3_active_genes_results.mat r lls params

[r_mef2, lls_mef2, params_mef2] = drosMLRank(drosexp, drosTF, 'mef2', drosexp.probes(active_genes));
r = r_mef2; lls = lls_mef2; params = params_mef2;
save results/dros_mlode_mef2_list_expro3_active_genes_results.mat r lls params

oderank.twi = drosRemoveDuplicateGenes(drosexp, r_twi);
oderank.mef2 = drosRemoveDuplicateGenes(drosexp, r_mef2);
