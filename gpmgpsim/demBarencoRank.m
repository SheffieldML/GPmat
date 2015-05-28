% DEMBARENCORANK Do ranking experiments on data from Barenco et al in Genome Biology.

% SHEFFIELDML

colordef white
[y, yvar, gene, times, scale, rawExp, rawVar] = gpsimLoadBarencoData;

load demBarenco1.mat

% Get the default options structure.
options = gpsimOptions;
options.includeNoise = 1;
% Fix one decay (from the fourth gene --- p21) to 0.8 hr^-1, and
% the corresponding sensitivity (see just after eqn 2 in the
% mathematical methods of Barenco et al.)
options.fix(1).index = 8;
options.fix(1).value = expTransform(0.8, 'xtoa');
options.fix(2).index = 9;
options.fix(2).value = expTransform(1, 'xtoa');


% initialise the model.
model.type = 'cgpsim';  % This new model type is a hack to run
                        % the model in a hierarchical manner.
                        % need to do this more elegantly later.
for i =1:3
  model.comp{i} = gpsimCreate(5, 1, times, y{i}, yvar{i}, options);
end

% Learn the model.
model = modelOptimise(model, [], [], 1, 3000);

% Each component of the model gives us a prediction given one set
% of replicates. The replicates are from an independent cell line,
% so the probabilistic assumption here is that the samples of f(t)
% are independent for each prediction (but inverse width, decay,
% basal rate and sensitivity parameters are all shared).
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
  
  % Info from Martino Paper:
  % 'True f' from Figure 3.
  % Got these figures with a ruler ...
  % Don't actually plot these below, but they are stored for reference
  truef = [0 1.6 2.6 2.5 2.6 1.6 0.9];
  truef = truef/sqrt(var(truef))*scalePred;
  
  
  % Figure 2(a) histograms;
  B = [2.6 1.5 0.5 0.2 1.35]; % From Martino paper ... but don't know the scale
  B = B/mean(B)*mean(model.comp{1}.B); % do a rough rescaling so
                                       % that the scales match.
  S = [3 0.8 0.7 1.8 0.7]/1.8; % From Martino paper ... but here we
                               % know the scale, because p21 is
                               % fixed to 1.
  D = [1.2 1.6 1.75 3.2 2.3]*0.8/3.2; % From Martino paper, again
                                      % we know the scale because
                                      % p21 is fixed to 0.8.
  
  % Martino f from Figure 2(b), again measured with a ruler.
  barencof = [0 2.7 3.9 2.3 1.5 1.6 1.4]/(1.8*mean(S))*mean(model.comp{1}.S);
  barencof = barencof/sqrt(var(barencof))*scalePred;
  
  
  figure, lin = plot(predt, predF, '-');
  hold on,
  bh = plot(predt, predF + 2*sqrt(varF), '--');
  bh = [bh plot(predt, predF - 2*sqrt(varF), '--')];
  lin = [lin plot(0:2:12, barencof, 'rx')];
  set(bh, 'lineWidth', 3);
  set(lin, 'lineWidth', 4);
  set(lin, 'markersize', 20);
  set(gca, 'fontname', 'arial', 'fontsize', 24, 'xlim', xlim)
  fileName = ['demBarenco1_profile' num2str(j)];
  print('-deps', ['../tex/diagrams/' fileName]);
  pos = get(gcf, 'paperposition')
  origpos = pos;
  pos(3) = pos(3)/2;
  pos(4) = pos(4)/2;
  set(gcf, 'paperposition', pos);
  lineWidth = get(gca, 'lineWidth');
  set(gca, 'lineWidth', lineWidth*2);
  %print('-dpng', ['../html/' fileName])
  set(gca, 'lineWidth', lineWidth);
  set(gcf, 'paperposition', origpos)

  
end



order = [1 5 3 4 2];

counter = 0;

% Plot first basal transcription rates.
figure
bar([model.comp{1}.B(order); B]', 0.6); colormap([0 0 0; 1 1 1]);
set(gca, 'xticklabel', {'DDB2', 'hPA26', 'TNFRSF20b', 'p21', 'BIK'})
fileName = ['demBarenco1_basal'];
print('-deps', ['../tex/diagrams/' fileName]);
pos = get(gcf, 'paperposition')
origpos = pos;
pos(3) = pos(3)/2;
pos(4) = pos(4)/2;
set(gcf, 'paperposition', pos);
lineWidth = get(gca, 'lineWidth');
set(gca, 'lineWidth', lineWidth*2);
print('-dpng', ['../html/' fileName])
set(gcf, 'paperposition', origpos)
set(gca, 'lineWidth', lineWidth);

% Plot the sensitivities.
figure
bar([model.comp{1}.S(order); S]', 0.6); colormap([0 0 0; 1 1 1]);
set(gca, 'xticklabel', {'DDB2', 'hPA26', 'TNFRSF20b', 'p21', ...
                    'BIK'})
fileName = ['demBarenco1_sensitivity'];
print('-deps', ['../tex/diagrams/' fileName]);
pos = get(gcf, 'paperposition')
origpos = pos;
pos(3) = pos(3)/2;
pos(4) = pos(4)/2;
set(gcf, 'paperposition', pos);
lineWidth = get(gca, 'lineWidth');
set(gca, 'lineWidth', lineWidth*2);
print('-dpng', ['../html/' fileName])
set(gcf, 'paperposition', origpos)
set(gca, 'lineWidth', lineWidth);

% Finally plot degradation rates.
figure
bar([model.comp{1}.D(order); D]', 0.6); colormap([0 0 0; 1 1 1]);
set(gca, 'xticklabel', {'DDB2', 'hPA26', 'TNFRSF20b', 'p21', ...
                    'BIK'})
fileName = ['demBarenco1_decay'];
print('-deps', ['../tex/diagrams/' fileName]);
pos = get(gcf, 'paperposition')
origpos = pos;
pos(3) = pos(3)/2;
pos(4) = pos(4)/2;
set(gcf, 'paperposition', pos);
lineWidth = get(gca, 'lineWidth');
set(gca, 'lineWidth', lineWidth*2);
print('-dpng', ['../html/' fileName])
set(gcf, 'paperposition', origpos)
set(gca, 'lineWidth', lineWidth);

