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

fprintf('\\begin{table}[p]\n');
fprintf('  \\centering\n');
drosPrintBootstrapMatrices(matrices{1}, 1:4)
fprintf('  \\caption{Bootstrap results: Twist global ChIP-chip}\n');
fprintf('  \\label{tab:bootstrap1}\n');
fprintf('\\end{table}\n\n');
fprintf('\\begin{table}[p]\n');
fprintf('  \\centering\n');
drosPrintBootstrapMatrices(matrices{2}, 1:4, 4)
fprintf('  \\caption{Bootstrap results: Twist global knock-outs}\n');
fprintf('  \\label{tab:bootstrap1}\n');
fprintf('\\end{table}\n\n');
fprintf('\\begin{table}[p]\n');
fprintf('  \\centering\n');
drosPrintBootstrapMatrices(matrices{3}, 1:4)
fprintf('  \\caption{Bootstrap results: Twist focused ChIP-chip}\n');
fprintf('  \\label{tab:bootstrap1}\n');
fprintf('\\end{table}\n\n');
fprintf('\\begin{table}[p]\n');
fprintf('  \\centering\n');
drosPrintBootstrapMatrices(matrices{4}, 1:4, 4)
fprintf('  \\caption{Bootstrap results: Twist focused knock-outs}\n');
fprintf('  \\label{tab:bootstrap1}\n');
fprintf('\\end{table}\n\n');
fprintf('\\begin{table}[p]\n');
fprintf('  \\centering\n');
drosPrintBootstrapMatrices(matrices{5}, 1:4)
fprintf('  \\caption{Bootstrap results: Mef2 global ChIP-chip}\n');
fprintf('  \\label{tab:bootstrap1}\n');
fprintf('\\end{table}\n\n');
fprintf('\\begin{table}[p]\n');
fprintf('  \\centering\n');
drosPrintBootstrapMatrices(matrices{6}, 1:4, 4)
fprintf('  \\caption{Bootstrap results: Mef2 global knock-outs}\n');
fprintf('  \\label{tab:bootstrap1}\n');
fprintf('\\end{table}\n\n');
fprintf('\\begin{table}[p]\n');
fprintf('  \\centering\n');
drosPrintBootstrapMatrices(matrices{7}, 1:4)
fprintf('  \\caption{Bootstrap results: Mef2 focused ChIP-chip}\n');
fprintf('  \\label{tab:bootstrap1}\n');
fprintf('\\end{table}\n\n');
fprintf('\\begin{table}[p]\n');
fprintf('  \\centering\n');
drosPrintBootstrapMatrices(matrices{8}, 1:4, 4)
fprintf('  \\caption{Bootstrap results: Mef2 focused knock-outs}\n');
fprintf('  \\label{tab:bootstrap1}\n');
fprintf('\\end{table}\n\n');
