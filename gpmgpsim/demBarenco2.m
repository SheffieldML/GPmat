% DEMBARENCO2 Run experiments on data from Barenco et al in Genome Biology.
% The raw data is pre-processed by MAS5.

% SHEFFIELDML
clear,clc
colordef white
[y, yvar, gene, times, scale] = gpsimLoadBarencoMASData;
iniTime = 0;
times = times + iniTime;

% Get the default options structure.
options = gpsimOptions;
options.includeNoise = 1;
% Prior information of the TF
options.proteinPrior = [0]';           % Assuming that f=0 at t=0 and t=12
options.proteinPriorTimes = [0]';
% Fix one decay (from the fourth gene --- p21) to 0.8 hr^-1, and
% the corresponding sensitivity (see just after eqn 2 in the
% mathematical methods of Barenco et al.)
if isfield(options, 'proteinPrior') && ~isempty(options.proteinPrior)
  options.fix(1).index = 9;
  options.fix(1).value = expTransform(0.8, 'xtoa');
  options.fix(2).index = 10;
  options.fix(2).value = expTransform(1, 'xtoa');
  options.fix(3).index = 2;              % RBF variance
  options.fix(3).value = expTransform(1, 'xtoa');
else
  options.fix(1).index = 8;
  options.fix(1).value = expTransform(0.8, 'xtoa');
  options.fix(2).index = 9;
  options.fix(2).value = expTransform(1, 'xtoa');
end


saveFigures = 0;

% initialise the model.
model.type = 'cgpsim';  % This new model type is a hack to run
                        % the model in a hierarchical manner.
                        % need to do this more elegantly later.
for i =1:3          %% 3 original
  model.comp{i} = gpsimCreate(5, 1, times, y{i}, yvar{i}, options);
end

% Learn the model.
model = modelOptimise(model, [], [], 1, 3000);

% Each component of the model gives us a prediction given one set
% of replicates. The replicates are from an independent cell line,
% so the probabilistic assumption here is that the samples of f(t)
% are independent for each prediction (but inverse width, decay,
% basal rate and sensitivity parameters are all shared).

% rescale the model parameters
modelB = model.comp{1}.B.*scale;
scaleModelB = modelB/mean(modelB);
modelS = model.comp{1}.S.*scale;
modelS = modelS/modelS(4);
order = [1 5 3 4 2];
geneNames = {'DDB2', 'hPA26', 'TNFRSF20b', 'p21', 'BIK'};
  
for j = 1:length(model.comp)
  predt = [linspace(0, 12, 100) 0:2:12]';
  predt = iniTime + predt;  

  ymean = reshape(ones(length(times),1)*y{j}(1,:), length(model.comp{j}.y), ...
                  1);
  
  % Generate predictions of the functions.
  % to do this we need to compute the K_xf portions of the kernel
  % (simXrbfKernCompute does this for us).

  if isfield(model.comp{j}, 'proteinPrior') && ~isempty(model.comp{j}.proteinPrior)
    if model.comp{j}.includeNoise
      for i=1:model.comp{j}.kern.comp{1}.numBlocks
        predTimeCell{i} = predt;
      end          
            
      Kxx = multiKernCompute(model.comp{j}.kern.comp{1}, predTimeCell, ...
                             model.comp{j}.timesCell);
      diagKxx = kernDiagCompute(model.comp{j}.kern.comp{1}, predTimeCell);     
      x = [model.comp{j}.proteinPrior; (model.comp{j}.y-ymean)];
      
      ind = 1:length(predTimeCell{1});
      for indBlock=1:model.comp{j}.kern.comp{1}.numBlocks
        K = Kxx(ind, :);
        diagK = diagKxx(ind,:);
        predFull{indBlock} = real(K*model.comp{j}.invK*x);
        varFull{indBlock} = real(diagK - sum(K'.*(model.comp{j}.invK*K'), 1)');
        ind = ind + length(predTimeCell{indBlock});
      end
    else
      for i=1:model.comp{j}.kern.numBlocks
        predTimeCell{i} = predt;
      end      
      
      Kxx = multiKernCompute(model.comp{j}.kern, predTimeCell, ...
                             model.comp{j}.timesCell);
      diagKxx = kernDiagCompute(model.comp{j}.kern, predTimeCell);
      x = [model.comp{j}.proteinPrior; (model.comp{j}.y-ymean)];
      
      ind = 1:length(predTimeCell{1});
      
      for indBlock=1:model.comp{j}.kern.numBlocks
        K = Kxx(ind, :);
        diagK = diagKxx(ind,:);
        predFull{indBlock} = real(K*model.comp{j}.invK*x);
        varFull{indBlock} = real(diagK - sum(K'.*(model.comp{j}.invK*K'), 1)');
        ind = ind + length(predTimeCell{indBlock});
      end
    end
    
    meanPredX = ones(length(predt),1)*(model.comp{j}.B./ ...
                      model.comp{j}.D);     
    predF = predFull{1};
    varF = varFull{1};
    predExprs = [];
    varExprs = [];
    for i = 1:model.comp{j}.numGenes
      predExprs(:,i) = meanPredX(:,i) + predFull{i+1};
      varExprs(:,i) = varFull{i+1};
    end
    predExprs(end-6:end,:) = [];   
    varExprs(end-6:end,:) = [];  
        
  else
    proteinKern = kernCreate(model.comp{1}.t, 'rbf'); 
    if model.comp{j}.includeNoise
      proteinKern.inverseWidth = ...
          model.comp{j}.kern.comp{1}.comp{1}.inverseWidth;
    else
      proteinKern.inverseWidth = ...
          model.comp{j}.kern.comp{1}.inverseWidth;
    end
    K = [];
  
    for i=1:model.comp{j}.kern.numBlocks
      if model.comp{j}.includeNoise
        K = [K; simXrbfKernCompute(model.comp{j}.kern.comp{1}.comp{i}, proteinKern, ...
                               model.comp{j}.t, predt)];
      else
        K = [K; simXrbfKernCompute(model.comp{j}.kern.comp{i}, proteinKern, ...
                               model.comp{j}.t, predt)];
      end
    end
  
    predF = real(K'*model.comp{j}.invK*(model.comp{j}.y-ymean));
    varF = real(kernDiagCompute(proteinKern, predt) - sum(K.*(model.comp{j}.invK*K), ...
                                                   1)');

    meanPredX = reshape(ones(length(predt),1)*(model.comp{j}.B./ ...
                      model.comp{j}.D), length(predt)*model.comp{j}.numGenes, ...
                        1);  
    
    if model.comp{j}.includeNoise
      Kxx = multiKernCompute(model.comp{j}.kern.comp{1}, predt, times);
      predX = meanPredX + real(Kxx*model.comp{j}.invK*(model.comp{j}.y-ymean));
      varX = real(kernDiagCompute(model.comp{j}.kern.comp{1}, predt) - sum(Kxx'.* ...
                                                    (model.comp{j}.invK*Kxx'), ...
                                                    1)');
    else
      Kxx = multiKernCompute(model.comp{j}.kern, predt, times);
      predX = meanPredX + real(Kxx*model.comp{j}.invK*(model.comp{j}.y-ymean));
      varX = real(kernDiagCompute(model.comp{j}.kern, predt) - sum(Kxx'.* ...
                                                    (model.comp{j}.invK*Kxx'), ...
                                                    1)');
    end
    numGenes = model.comp{j}.numGenes;
    numTimes = length(predX)/numGenes;
    predExprs = reshape(predX, numTimes, numGenes);
    predExprs(end-6:end,:) = [];
    varExprs = reshape(varX, numTimes, numGenes);
    varExprs(end-6:end,:) = [];  
  end
  
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
  B = B/mean(B)*mean(modelB); % do a rough rescaling so
                                       % that the scales match.
  S = [3 0.8 0.7 1.8 0.7]/1.8; % From Martino paper ... but here we
                               % know the scale, because p21 is
                               % fixed to 1.
  D = [1.2 1.6 1.75 3.2 2.3]*0.8/3.2; % From Martino paper, again
                                      % we know the scale because
                                      % p21 is fixed to 0.8.
  
  % Martino f from Figure 2(b), again measured with a ruler.
  
  barencof = [0.0000000 200.5201100 355.5216125 205.7574913 135.0911372 ...
              145.1080997 130.7046969; 0.0000000 184.0994134 308.4759200 ...
              232.1775328 153.6595161 85.7272235  168.0910562; 0.0000000 ...
              230.2262511 337.5994811 276.9416540 164.5044287 127.8653452 173.6112139];
  
%  barencof = barencof/(1.8*mean(S))*mean(modelS);
  barencof = barencof./(sqrt(var(barencof,0,2))*ones(1,7))*scalePred;
  
  figure,
  lin = plot(predt, predF, '-');
  hold on,
  bh = plot(predt, predF + 2*sqrt(varF), '--');
  bh = [bh plot(predt, predF - 2*sqrt(varF), '--')];
%  lin = [lin plot(times, barencof(j,:), 'rx')];
  titleText = 'Inferred p53 protein';
  title(titleText,'fontsize',20);
  set(bh, 'lineWidth', 3);
  set(lin, 'lineWidth', 4);
  set(lin, 'markersize', 20);
  set(gca, 'fontname', 'arial', 'fontsize', 24, 'xlim', xlim)
  if saveFigures==1
    fileName = ['demBarenco1_profile' num2str(j)];
    print('-deps', ['./results/' fileName]);
    pos = get(gcf, 'paperposition');
    origpos = pos;
    pos(3) = pos(3)/2;
    pos(4) = pos(4)/2;
    set(gcf, 'paperposition', pos);
    lineWidth = get(gca, 'lineWidth');
    set(gca, 'lineWidth', lineWidth);
    print('-dpng', ['./results/' fileName])
    set(gca, 'lineWidth', lineWidth);
    set(gcf, 'paperposition', origpos);
  end
  
  for index = 1:model.comp{j}.numGenes
    figure;
%     subplot(3,2,index);
    lin = plot(predt, predExprs(:,index), '-');
    hold on,
    bh = plot(predt, predExprs(:,index)+2*real(sqrt(varExprs(:,index))), '--');
    bh = [bh plot(predt, predExprs(:,index)-2*real(sqrt(varExprs(:,index))), '--')];
    lin1 = errorbar(times, y{j}(:,index), 2*sqrt(yvar{j}(:,index)), ...
                    'rx');
    lin = [lin plot(times, y{j}(:,index), 'rx')];
    titleText = ['gene ' geneNames{order(index)} ' mRNA'];
    title(titleText,'fontsize', 20);
    textB = ['B = ' num2str(scaleModelB(index))];
    textD = ['D = ' num2str(model.comp{j}.D(index))];
    textS = ['S = ' num2str(modelS(index))];
    texPos = ylim;
    scalePos = (texPos(2) - texPos(1))/14;
    texh = text(9, texPos(2)-scalePos, textB);
    texh = [texh text(9, texPos(2)-2*scalePos,textD)];
    texh = [texh text(9, texPos(2)-3*scalePos, textS)];
    set(bh, 'lineWidth', 3);
    set(lin, 'lineWidth', 4);
    set(lin, 'markersize', 20);
    set(lin1, 'lineWidth', 2);
    set(texh, 'fontsize', 16);
    set(gca, 'fontname', 'arial', 'fontsize', 24, 'xlim', [min(predt) max(predt)])
  
    if saveFigures==1
      fileName = ['demBarenco1_ExprsProfile_Rep' num2str(j) '_Gene' num2str(index)];
      print('-deps', ['./results/' fileName]);
      pos = get(gcf, 'paperposition');
      origpos = pos;
      pos(3) = pos(3);
      pos(4) = pos(4);
      set(gcf, 'paperposition', pos);
      lineWidth = get(gca, 'lineWidth');
      set(gca, 'lineWidth', lineWidth*2);
      print('-dpng', ['./results/' fileName]);
      set(gca, 'lineWidth', lineWidth);
      set(gcf, 'paperposition', origpos);
    end
  end
     
end

counter = 0;

% Plot first basal transcription rates.
figure
bar([modelB(order); B]', 0.6); colormap([0 0 0; 1 1 1]);
set(gca, 'xticklabel', {'DDB2', 'hPA26', 'TNFRSF20b', 'p21', 'BIK'})
if saveFigures==1
  fileName = ['demBarenco1_basal'];
  print('-deps', ['./results/' fileName]);
  pos = get(gcf, 'paperposition');
  origpos = pos;
  pos(3) = pos(3)/2;
  pos(4) = pos(4)/2;
  set(gcf, 'paperposition', pos);
  lineWidth = get(gca, 'lineWidth');
  set(gca, 'lineWidth', lineWidth*2);
  print('-dpng', ['./results/' fileName])
  set(gcf, 'paperposition', origpos)
  set(gca, 'lineWidth', lineWidth);
end

% Plot the sensitivities.
figure
bar([modelS(order); S]', 0.6); colormap([0 0 0; 1 1 1]);
set(gca, 'xticklabel', {'DDB2', 'hPA26', 'TNFRSF20b', 'p21', ...
                    'BIK'})
if saveFigures==1
  fileName = ['demBarenco1_sensitivity'];
  print('-deps', ['./results/' fileName]);
  pos = get(gcf, 'paperposition');
  origpos = pos;
  pos(3) = pos(3)/2;
  pos(4) = pos(4)/2;
  set(gcf, 'paperposition', pos);
  lineWidth = get(gca, 'lineWidth');
  set(gca, 'lineWidth', lineWidth*2);
  print('-dpng', ['./results/' fileName])
  set(gcf, 'paperposition', origpos)
  set(gca, 'lineWidth', lineWidth);
end

% Finally plot degradation rates.
figure
bar([model.comp{1}.D(order); D]', 0.6); colormap([0 0 0; 1 1 1]);
set(gca, 'xticklabel', {'DDB2', 'hPA26', 'TNFRSF20b', 'p21', ...
                    'BIK'})
if saveFigures==1
  fileName = ['demBarenco1_decay'];
  print('-deps', ['./results/' fileName]);
  pos = get(gcf, 'paperposition');
  origpos = pos;
  pos(3) = pos(3)/2;
  pos(4) = pos(4)/2;
  set(gcf, 'paperposition', pos);
  lineWidth = get(gca, 'lineWidth');
  set(gca, 'lineWidth', lineWidth*2);
  print('-dpng', ['./results/' fileName])
  set(gcf, 'paperposition', origpos)
  set(gca, 'lineWidth', lineWidth);
end
