% DEMRUNRANKINGS Runs the rankings presented in the paper
% FORMAT
% DESC Runs the rankings presented in the paper (may take a really long time...)
%
% COPYRIGHT : Antti Honkela, 2009

% SHEFFIELDML

fprintf('Warning: running this demo may take *SEVERAL DAYS OR WEEKS*\n');
fprintf('even on a very fast computer if the rankings need to be rerun.\n')
fprintf('Press Ctrl-C to cancel or any other key to continue...\n')
pause;

drosLoadData;
z_scores = drosexp.fitmean ./ sqrt(drosexp.fitvar);
active_genes = (mean(z_scores, 2) >= 1.8);

tfnames = {'twi', 'mef2'};

% Load rankings if the files exist, else run the rankings (and save)
for k=1:length(tfnames),
  tf = tfnames{k};
  listfile = 'expro3_active_genes';
  stfname = sprintf('results/dros_gpdisim_%s_list_%s_results.mat', tf, listfile);
  if exist(stfname, 'file'),
    stres.(tf) = load(stfname);
  else
    stres.(tf) = drosScoreTFTargetList(tf, drosexp.probes(active_genes), []);
    save(stfname, '-struct', stres.(tf));
  end

  if strcmp(tf, 'twi'),
    ssfname = sprintf('results/dros_gpdisim_subsample_%s_list_%s_results.mat', tf, listfile);
    if exist(ssfname, 'file'),
      ssres.(tf) = load(ssfname);
    else
      ssres.(tf) = drosScoreTFTargetList(tf, drosexp.probes(active_genes), [], 1, 0, [2:6, 8, 10]);
    end
  end
  
  basefile = 'basetargets_kosher';
  mtfname = sprintf('results/dros_gpdisimmt_%s_list_%s_dros_%s_%s_results.mat', tf, listfile, tf, basefile);
  if exist(mtfname, 'file'),
    mtres.(tf) = load(mtfname);
  else
    [dummy, K] = sort(stres.(tf).ll, 'descend');
    ranking = stres.(tf).targets(K);
    % Exclude the TF itself from the ranking
    I = ~strcmp(drosTF.probes.(tf), ranking);
    ranking = ranking(I);
  
    mtres.(tf) = drosScoreTFTargetListMT(tf, drosexp.probes(active_genes), ...
					 ranking(1:5), []);
    save(mtfname, '-struct', mtres.(tf));
  end
end

% Process the results
drosexp2 = drosexp;
drosexp2.pctiles = drosexp2.pctiles(active_genes, :, :);
drosexp2.genes = drosexp2.genes(active_genes, :);
drosexp2.probes = drosexp2.probes(active_genes, :);
drosexp2.symbols = drosexp2.symbols(active_genes, :);
drosexp2.mean = drosexp2.mean(active_genes, :);
drosexp2.se = drosexp2.se(active_genes, :);
drosexp2.fitmean = drosexp2.fitmean(active_genes, :);
drosexp2.fitvar = drosexp2.fitvar(active_genes, :);

for k=1:length(tfnames),
  tf = tfnames{k};
  AG = drosFindGeneinds(drosexp, stres.(tf).targets, 0, 1);
  [ll_sort, K] = sort(stres.(tf).ll(active_genes(AG)), 'descend');
  indrank.(tf) = drosRemoveDuplicateGenes(drosexp, AG(K));

  if k==1,
    inds0 = [2:6 8 10];
    ss_inds = [inds0, inds0+12, inds0+24];
    z_scores2 = drosexp.fitmean(:, ss_inds) ./ sqrt(drosexp.fitvar(:, ss_inds));
    active_genes2 = (mean(z_scores2, 2) >= 1.8);

    AG = drosFindGeneinds(drosexp, ssres.(tf).targets, 0, 1);
    [actives, I1, J1] = intersect(AG, find(active_genes2));
    [ll_sort, K] = sort(ssres.(tf).ll(I1), 'descend');
    ssrank.(tf) = drosRemoveDuplicateGenes(drosexp, actives(K));
  end

  AG = drosFindGeneinds(drosexp, mtres.(tf).targets, 0, 1);
  [ll_sort, K] = sort(mtres.(tf).ll(active_genes(AG)), 'descend');
  disimrank.(tf) = drosRemoveDuplicateGenes(drosexp, AG(K));

  CR = drosCorrelationRank(drosexp, drosTF, tf);
  corrrank.(tf) = drosRemoveDuplicateGenes(drosexp, CR(active_genes(CR)));

  if exist('drosmutant'),
    if isfield(drosmutant, tf),
      J = drosFindGeneinds(drosexp, drosmutant.(tf).genes, 1);
      mutarank.(tf) = drosRemoveDuplicateGenes(drosexp, J(J ~= 0));
    else
      mutarank.(tf) = [];
    end
  end
end
