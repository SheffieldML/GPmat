% DROSDOBOOTSTRAP Perform bootstrap sampling to assess significance of differences in rankings
% FORMAT
% DESC Perform bootstrap sampling to assess significance of differences in rankings
%
% COPYRIGHT : Antti Honkela, 2009

% DISIMRANK

demRunRankings

thresholds = [5000];
REPEATS = 100000;
tt = [20, 100, 250];
matrices = {};
data = {};

tfs = {'twi', 'mef2'};

for i=1:2,
  tf = tfs{i};
  [a, b] = drosBootstrapEvaluation(drosexp, drosinsitu, {indrank.(tf), disimrank.(tf), corrrank.(tf), mutarank.(tf)}, chip_validation.(tf), REPEATS, tt, 0);
  a, mean(b)
  matrices = [matrices, {a}];
  data = [data, {b}];

  [a, b] = drosBootstrapEvaluation(drosexp, drosinsitu, {indrank.(tf), disimrank.(tf), corrrank.(tf), mutarank.(tf)}, mutant_validation.(tf), REPEATS, tt, 0);
  a, mean(b)
  matrices = [matrices, {a}];
  data = [data, {b}];

  [a, b] = drosBootstrapEvaluation(drosexp, drosinsitu, {indrank.(tf), disimrank.(tf), corrrank.(tf), mutarank.(tf)}, chip_validation.(tf), REPEATS, tt, 1);
  a, mean(b)
  matrices = [matrices, {a}];
  data = [data, {b}];

  [a, b] = drosBootstrapEvaluation(drosexp, drosinsitu, {indrank.(tf), disimrank.(tf), corrrank.(tf), mutarank.(tf)}, mutant_validation.(tf), REPEATS, tt, 1);
  a, mean(b)
  matrices = [matrices, {a}];
  data = [data, {b}];
end
%i = 2; [a, b] = drosBootstrapEvaluation(drosexp2, drosinsitu, {indids{i}, disimids{i}, corrids{i}, mutaids{i}}, validation.(tfnames{i}), 1000, [20 100 250], 1); a, mean(b)
%i = 2; [a, b] = drosBootstrapEvaluation(drosexp2, drosinsitu, {indids{i}, disimids{i}, corrids{i}, mutaids{i}}, validation.(tfnames{i}), 1000, [20 100 250], 0); a, mean(b)

save /share/work/ahonkela/bootstrap_results2.mat matrices data

fprintf('\\multicolumn{15}{c}{Twist global ChIP-chip} \\\\\n')
drosPrintBootstrapMatrices(matrices{1}, 1:4)
fprintf('\\multicolumn{15}{c}{Twist global knock-outs} \\\\\n')
drosPrintBootstrapMatrices(matrices{2}, 1:4, 4)
fprintf('\\multicolumn{15}{c}{Twist focused ChIP-chip} \\\\\n')
drosPrintBootstrapMatrices(matrices{3}, 1:4)
fprintf('\\multicolumn{15}{c}{Twist focused knock-outs} \\\\\n')
drosPrintBootstrapMatrices(matrices{4}, 1:4, 4)
fprintf('\\multicolumn{15}{c}{Mef2 global ChIP-chip} \\\\\n')
drosPrintBootstrapMatrices(matrices{5}, 1:4)
fprintf('\\multicolumn{15}{c}{Mef2 global knock-outs} \\\\\n')
drosPrintBootstrapMatrices(matrices{6}, 1:4, 4)
fprintf('\\multicolumn{15}{c}{Mef2 focused ChIP-chip} \\\\\n')
drosPrintBootstrapMatrices(matrices{7}, 1:4)
fprintf('\\multicolumn{15}{c}{Mef2 focused knock-outs} \\\\\n')
drosPrintBootstrapMatrices(matrices{8}, 1:4, 4)
