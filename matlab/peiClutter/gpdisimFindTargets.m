% function model = gpdisimDemDrosoph(drosexp, drosTF, tf),
clear, clc; close all;
expNo = 1;
type = 'DrosSingleGene';
drosLoadMef2Data;
targetsFull = drosFindTargets(drosmef2chip);
warning('off');
saveFigures = 0;

tf = 'mef2';
% targets = drosTargets(tf); 
tflabel = drosTF.labels(strcmp(tf, drosTF.names));

tfName = tf;
tfName(1) = upper(tfName(1));
fileName = ['dem' tfName type num2str(expNo)];
nTestGenes = 50;
nSelection = 20;
ll = -1000*ones(1,nTestGenes);

for targetInd = 1:nTestGenes
% targets = ([targetsFull(1:6); targetsFull(8:10)])';
targets = targetsFull(targetInd);
targetGene = drosGetGeneinds(drosexp, targets);
if ~isempty(targetGene)
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

  try,
% Learn the model.
    model = modelOptimise(model, [], [], 1, 3000);
    modelFull{targetInd} = model;
    ll(targetInd) = cgpdisimLogLikelihood(model);
  catch
    fprintf('%d-th gene failed on the model optimisation!', targetInd);
  end
end
end

[selectLl selectInd] = sort(ll, 'descend');

for targetInd = 11:nSelection
  model = modelFull{selectInd(targetInd)};

% Plot
% figure(1); clf; plot_tf_exps(exp_struct, expse_struct);
% figure(2); clf; plot_expros(genes, exp_struct, expse_struct, genenames);
% figure(3); clf;

  numGenes = model.comp{1}.numGenes;

  for j = 1:length(model.comp)

  % Generate predictions of the functions.
  % to do this we need to compute the K_xf portions of the kernel
  % (simXrbfKernCompute does this for us).
    predt = [linspace(1, 12, 100) model.comp{j}.t']';
    proteinKern = kernCreate(model.comp{j}.t, 'sim');
    proteinKern.inverseWidth = model.comp{j}.kern.comp{1}.inverseWidth;
    proteinKern.decay = model.comp{j}.delta;
    proteinKern.variance = model.comp{j}.kern.comp{2}.di_variance;  
    K = simXrbfKernCompute(proteinKern, model.comp{j}.kern.comp{1}, ...
                               predt, model.comp{j}.t);
    for i=2:model.comp{j}.kern.numBlocks
      blockK =  disimXsimKernCompute(model.comp{j}.kern.comp{i}, proteinKern, ...
                               model.comp{j}.t, predt);
      K = [K blockK'];    
    end
    predF = K*model.comp{j}.invK*model.comp{j}.y;
    varF = abs( kernDiagCompute(proteinKern, predt) - sum(K'.* ...
                                                   (model.comp{j}.invK*K'))'); ...
           %???
    
  % Predicted Gene Expressions
    Kxx = multiKernCompute(model.comp{j}.kern, predt, model.comp{j}.t);
    predX = real(Kxx*model.comp{j}.invK*model.comp{j}.y);
    varX = real(kernDiagCompute(model.comp{j}.kern, predt) - sum(Kxx'.* ...
                           (model.comp{j}.invK*Kxx'), 1)');
    
  % Take out predictions at data points.
  % Use them to get the scale for the other data.
    numData = length(predX)/(numGenes+1);
    predExprs = reshape(predX, numData, numGenes+1);
    meanExprs = ones(numData, 1)*mean(predExprs);
    scaleExprs = ones(numData, 1)*sqrt(var(predExprs));
    predExprs = meanExprs+(predExprs - meanExprs)./scaleExprs;
    predExprs(end-length(model.comp{j}.t):end,:) = [];
    varExprs = reshape(varX, numData, numGenes+1);
    varExprs = varExprs./scaleExprs./scaleExprs;
    varExprs(end-length(model.comp{j}.t):end,:) = []; 
    
    numTimes = length(model.comp{j}.t);
    scaleY{j} = ones(numTimes, 1)*mean(y{j})+(y{j}-ones(numTimes,1)*mean(y{j}))./(ones(numTimes,1)*sqrt(var(y{j})));

    dataF = predF(end-length(model.comp{j}.t):end);
    dataVarF = varF(end-length(model.comp{j}.t):end);
    predF(end-length(model.comp{j}.t):end) = [];
    varF(end-length(model.comp{j}.t):end) = [];
    predt(end-length(model.comp{j}.t):end) = [];
    scalePred = sqrt(var(dataF));
  
%     figure(targetInd);
%     subplot(length(model.comp), 1, j);
%     lin = plot(predt, predF, '-');
%     hold on,
%     bh = plot(predt, predF + 2*sqrt(varF), '--');
%     bh = [bh plot(predt, predF - 2*sqrt(varF), '--')];
%     hold off
%     ylabel(['Replica' num2str(j)]);
%     set(bh, 'lineWidth', 2);
%     set(lin, 'lineWidth', 3);
%   %set(lin, 'markersize', 20);
%     set(gca, 'fontname', 'arial', 'fontsize', 20, 'xlim', [min(model.comp{j}.t) ...
%                                                     max(model.comp{j} ...
%                                                     .t)]);
  
    for index = 1:numGenes+1
      figure(nTestGenes+targetInd);
      subplot(3, numGenes+1,index+(j-1)*(numGenes+1));
      lin = plot(predt, predExprs(:,index), '-');
      hold on,
      bh = plot(predt, predExprs(:,index)+2*real(sqrt(varExprs(:,index))), '--');
      bh = [bh plot(predt, predExprs(:,index)-2*real(sqrt(varExprs(:,index))), ...
                    '--')];
      lin = [lin plot(times, scaleY{j}(:,index), 'rx')];
      titleText = [num2str(index) '-th mRNA in Rep' num2str(j) ':Likelihood ' ...
                          'of' num2str(selectLl(targetInd))];
      title(titleText);
      set(bh, 'lineWidth', 1);
      set(lin, 'lineWidth', 2);
      set(lin, 'markersize', 5);
      set(gca, 'fontname', 'arial', 'fontsize', 24, 'xlim', [min(predt) max(predt)])
    end
  end
  if saveFigures==1
%       figure(targetInd)
%       saveFileName = [fileName '_TFProfile_Gene' num2str(targetInd)];
%       print('-deps', ['./results/' saveFileName]);
%       pos = get(gcf, 'paperposition');
%       origpos = pos;
%       pos(3) = pos(3);
%       pos(4) = pos(4);
%       set(gcf, 'paperposition', pos);
%       lineWidth = get(gca, 'lineWidth');
%       set(gca, 'lineWidth', lineWidth*2);
%       print('-dpng', ['./results/' saveFileName]);
%       set(gca, 'lineWidth', lineWidth);
%       set(gcf, 'paperposition', origpos);
      
      figure (nTestGenes+targetInd);
      saveFileName = [fileName 'predExprs_profile' num2str(targetInd)];
      print('-deps', ['./results/' saveFileName]);
      pos = get(gcf, 'paperposition');
      origpos = pos;
      pos(3) = pos(3)/2;
      pos(4) = pos(4)/2;
      set(gcf, 'paperposition', pos);
      lineWidth = get(gca, 'lineWidth');
      set(gca, 'lineWidth', lineWidth);
      print('-dpng', ['./results/' saveFileName])
      set(gca, 'lineWidth', lineWidth);
      set(gcf, 'paperposition', origpos);
  end
end

save(fileName);