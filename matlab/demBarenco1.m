% DEMBARENCO1 Run experiments on data from Barenco et al in Genome Biology.

% GPSIM

colordef white

if exist('./data/demBarenco1.mat') == 2
  load('./data/demBarenco1.mat');
else

  % These excel files include results processed directly from the
  % cel files using the mmgMOS algorithm (Xuejun's code).

  % These are the expression levels.
  [numeric1, txt1] = xlsread('./data/resultsMartino_exprs.xls');
  headTxt1 = txt1(1, 2:end);
  tagTxt1 = txt1(2:end, 1);

  % These are the standard deviations.
  [numeric2, txt2] = xlsread('./data/resultsMartino_se.xls');
  headTxt2 = txt2(1, 2:end);
  tagTxt2 = txt2(2:end, 1);
  
  if(any(~strcmp(tagTxt2(:), tagTxt1(:))))
    error('Two files are not in same order');
  end
  if(any(~strcmp(headTxt2(:), headTxt1(:))))
    error('Two files are not in same order');
  end
  
  clear gene, clear ind
  % Gene IDs
  % DDB2
  gene{1, 1} = '203409_at';
  gene{1, 2} = 'DDB2';
  % BIK
  gene{2, 1} = '205780_at';
  gene{2, 2} = 'BIK';
  % TNFRSF10b (other tags include 209294_x_at and 210405_x_at)
  gene{3, 1} = '209295_at';
  gene{3, 2} = 'TNFRSF10b';
  % p21 --- we think this is CIp1/p21
  gene{4, 1} = '202284_s_at';
  gene{4, 2} = 'CIp1/p21';
  % p26 --- named as sesn1 in the platform.
  gene{5, 1} = '218346_s_at';
  gene{5, 2} = 'p26 sesn1';
  
  for i = 1:length(gene)
    match = find([strcmp(gene{i, 1}, tagTxt1(:))]);
    if length(match)~=1
      error('Too many or too few matches.');
    else
      ind(i) = match(1);
    end
  end
  order = [1 4 5 6 7 2 3 8 11 12 13 14 9 10 15 18 19 20 21 16 17]; 

  
  % Perform some normalisation.
  % Make sure that the average for each slide in log space is the
  % same.
  mVal = zeros(size(mean(numeric1)));
  mVal = mVal - mean(mVal);
  rawExp = numeric1(ind, order)';
  for i = 1:size(rawExp, 2)
    rawExp(:, i) = rawExp(:, i) - mVal';
  end
  
  rawVar = numeric2(ind, order)';
  rawVar = rawVar.*rawVar; % convert standard deviations to variances.

  yFull = exp(rawExp + rawVar/2);  % Logs are normally distributed
                               % ... recover mean in exp space.
  yFullVar = (exp(rawVar)-1).*exp(2*rawExp + rawVar); % Logs are
                                                  % normally
                                                  % distributed
                                                  % ... recover
                                                  % variance in exp
                                                  % space.

  
  % Rescale so that average standard deviation of curves is 1.
  scale = mean(sqrt(var(yFull)));
  yFull = yFull/scale;
  yFullVar = yFullVar/(scale*scale);
  y{1} = yFull(1:7, :);
  y{2} = yFull(8:14, :);
  y{3} = yFull(15:21, :);
  yvar{1} = yFullVar(1:7, :);
  yvar{2} = yFullVar(8:14, :);
  yvar{3} = yFullVar(15:21, :);
  times = [0 2 4 6 8 10 12]';
  save('./data/demBarenco1.mat', 'y', 'yvar', 'gene', 'times', 'scale', 'rawVar', 'rawExp');
end

% Get the default options structure.
options = gpsimOptions;
options.includeNoise = 1;
% Fix one decay (from the fourth gene --- p21) to 0.8 hr^-1, and
% the corresponding sensitivity (see just after eqn 2 in the
% mathematical methods of Barenco et al.)
options.fix(1).index = 8;
options.fix(1).value = negLogLogitTransform(0.8, 'xtoa');;
options.fix(2).index = 9;
options.fix(2).value = negLogLogitTransform(1, 'xtoa');;


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
  plot(predt, predF + 2*sqrt(varF), '--')
  plot(predt, predF - 2*sqrt(varF), '--')
  lin = [lin plot(0:2:12, barencof, 'x')];
  set(lin, 'lineWidth', 2);
  set(lin, 'markersize', 10);
  fileName = ['demBarenco1_profile' num2str(j)];
  print('-deps', ['../tex/diagrams/' fileName]);
  pos = get(gcf, 'paperposition')
  origpos = pos;
  pos(3) = pos(3)/2;
  pos(4) = pos(4)/2;
  set(gcf, 'paperposition', pos);
  lineWidth = get(gca, 'lineWidth');
  set(gca, 'lineWidth', lineWidth*2);
  print('-dpng', ['../html/' fileName])
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

