function model = gpdisimDemDrosoph(drosexp, drosTF, tf),

targets = drosTargets(tf);
tflabel = drosTF.labels(strcmp(tf, drosTF.names));

genes = [tflabel, targets];

genenames = genes;

[y, yvar, gene, times, scale, rawExp, rawVar] = gpdisimGetDrosData(drosexp, genes);

% Get the default options structure.
options = gpsimOptions;
options.includeNoise = 1;
% Fix one decay (from the fourth gene --- p21) to 0.8 hr^-1, and
% the corresponding sensitivity (see just after eqn 2 in the
% mathematical methods of Barenco et al.)
options.fix(1).index = 2;
options.fix(1).value = expTransform(1, 'xtoa');
options.fix(2).index = 6;
options.fix(2).value = expTransform(1, 'xtoa');
%options.fix(2).index = 1;
%options.fix(2).value = expTransform(.3, 'xtoa');;


% initialise the model.
model.type = 'cgpdisim';  % This new model type is a hack to run
                        % the model in a hierarchical manner.
                        % need to do this more elegantly later.
for i =1:3          %% 3 original
  model.comp{i} = gpdisimCreate(length(targets), 1, times, y{i}, yvar{i}, options);
end

% Learn the model.
model = modelOptimise(model, [], [], 1, 3000);


return


%%% THIS NEEDS TO BE FIXED!


% Plot
figure(1); clf; plot_tf_exps(exp_struct, expse_struct);
figure(2); clf; plot_expros(genes, exp_struct, expse_struct, genenames);
figure(3); clf;

for j = 1:length(model.comp)

  % Generate predictions of the functions.
  % to do this we need to compute the K_xf portions of the kernel
  % (simXrbfKernCompute does this for us).
  predt = [linspace(-2, 14, 100) 0:2:12]';
  proteinKern = kernCreate(model.comp{1}.t, 'rbf');
  proteinKern.inverseWidth = ...
      model.comp{j}.kern.comp{1}.inverseWidth;
  K = [];
  for i=1:model.comp{j}.kern.numBlocks
    K = [K; simXrbfKernCompute(model.comp{j}.kern.comp{i}, proteinKern, ...
                               model.comp{j}.t, predt)];
    
  end
  predF = K'*model.comp{j}.invK*model.comp{j}.y;
  varF = kernDiagCompute(proteinKern, predt) - sum(K.*(model.comp{j}.invK*K), 1)';

  % Take out predictions at data points.
  % Use them to get the scale for the other data.
  dataF = predF(end-6:end);
  dataVarF = varF(end-6:end);
  predF(end-6:end) = [];
  varF(end-6:end) = [];
  predt(end-6:end) = [];
  scalePred = sqrt(var(dataF));

  subplot(length(model.comp), 1, j);
  lin = plot(predt, predF, '-');
  hold on,
  bh = plot(predt, predF + 2*sqrt(varF), '--');
  bh = [bh plot(predt, predF - 2*sqrt(varF), '--')];
  hold off
  %set(bh, 'lineWidth', 3);
  %set(lin, 'lineWidth', 4);
  %set(lin, 'markersize', 20);
  %set(gca, 'fontname', 'arial', 'fontsize', 24, 'xlim', xlim)
end
