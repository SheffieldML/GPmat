

genes1 = {'FBgn0003117', 'FBgn0029082', 'FBgn0033459', 'FBgn0036549', 'FBgn0051368'};

genenames = {'pannier', 'hibris', 'CG12744', 'CG10516', 'CG31368'};

genes = genes1;

%[y, yvar, gene, times, scale, rawExp, rawVar] = gpsimGetDrosData(exp_struct, expse_struct, genes);
[y, yvar, gene, times, scale, rawExp, rawVar] = gpsimGetDrosData(pcdata, [], genes);

% Get the default options structure.
options = gpsimOptions;
options.includeNoise = 1;
% Fix one decay (from the fourth gene --- p21) to 0.8 hr^-1, and
% the corresponding sensitivity (see just after eqn 2 in the
% mathematical methods of Barenco et al.)
options.fix(1).index = 3;
options.fix(1).value = expTransform(1, 'xtoa');;
%options.fix(2).index = 1;
%options.fix(2).value = expTransform(.3, 'xtoa');;


% initialise the model.
model.type = 'cgpsim';  % This new model type is a hack to run
                        % the model in a hierarchical manner.
                        % need to do this more elegantly later.
for i =1:3          %% 3 original
  model.comp{i} = gpsimCreate(length(genes), 1, times, y{i}, yvar{i}, options);
end

% Learn the model.
model = modelOptimise(model, [], [], 1, 3000);


% Plot
figure(1);
plot_tf_exps(exp_struct, expse_struct);
subplot(5, 1, 1);
title('TF profiles')
figure(2);
plot_expros(genes, exp_struct, expse_struct, genenames);
subplot(5, 1, 1);
title('Target profiles')
figure(3);

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
  set(gca, 'fontname', 'arial', 'fontsize', 24, 'xlim', xlim)
end
