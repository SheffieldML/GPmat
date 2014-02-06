% DEMGPDISIMMEF2 Run experiments on Mef2 data. The raw data is pre-processed by the PUMA package.

% SHEFFIELDML

%/~
% path(path, '/local/Matlab/underDevelopment/GPSIM016');
%~/

clear; close all;
expNo = 2;
type = 'Dros';
warning('off');
saveFigures = 0;

tf = 'mef2';

%/~
if exist('./data/mef2Data.mat') == 2 
%~/
load('./data/mef2Data.mat');
%/~
else
  
  drosLoadMef2Data;
  targetsFull = drosFindTargets(drosmef2chip);
  % targets = ([targetsFull(1:6); targetsFull(8:10)])';
  % selection = [8 34 19 26 2];
  selection = [25 34 19 37 21 40];
  targets = targetsFull(selection)';
  
  % targets = drosTargets(tf); 
  tflabel = drosTF.labels(strcmp(tf, drosTF.names));
  
  genes = [tflabel, targets];
  
  [y, yvar, gene, times, scale, rawExp, rawVar] = gpdisimGetDrosData(drosexp, ...
                                                    genes);
  
  save('./data/mef2Data.mat', 'y', 'yvar', 'gene', 'times', 'scale', 'rawVar', ...
       'rawExp', 'genes', 'targets');
  
end
%~/
genenames = genes;

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

tfName = tf;
tfName(1) = upper(tfName(1));
fileName = ['dem' tfName type num2str(expNo)];
save(fileName);

% Plot
% figure(1); clf; plot_tf_exps(exp_struct, expse_struct);
% figure(2); clf; plot_expros(genes, exp_struct, expse_struct, genenames);
% figure(3); clf;

numGenes = model.comp{1}.numGenes;
genenames{2} = 'Rya-r44F';
genenames{4} = 'ttk';

for j = 1:length(model.comp)

  % Generate predictions of the functions.
  % to do this we need to compute the K_xf portions of the kernel
  % (simXrbfKernCompute does this for us).
  predt = [1:0.1:12 model.comp{j}.t']';
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
  
  ymean = reshape(ones(length(times),1)*y{j}(1,:), length(model.comp{j}.y), ...
                                                     1);
%  ymean = mean(model.comp{j}.y);

  predF = K*model.comp{j}.invK*(model.comp{j}.y-ymean);
  varF = kernDiagCompute(proteinKern, predt) - sum(K'.* ...
                                                   (model.comp{j}.invK* ...
                                                    K'))';
%  varF = proteinKern.variance - sum(K'.*(model.comp{j}.invK*K'))';
  
%   % predict gene via the model
%   predF(end-length(model.comp{j}.t)+1:end) = [];
%   varF(end-length(model.comp{j}.t)+1:end) = [];
%   predt(end-length(model.comp{j}.t)+1:end) = [];  
%   model.comp{j}.mapt = predt;
%   model.comp{j}.g = predF;
%   model.comp{j}.numMapPts = length(model.comp{j}.mapt);
%   model.comp{j}.step = 0.1;
%    
%   model.comp{j}.times_index = [];
%   for i = 1:length(times)
%     model.comp{j}.times_index(i) = find((times(i) - model.comp{j}.mapt)==0);
%   end
%   
%   model.comp{j} = gpsimMapUpdateYpred(model.comp{j});
%   
%   predExprs = zeros(length(predt), model.comp{j}.numGenes+1);
%   varExprs = zeros(length(predt), model.comp{j}.numGenes+1);  
%   
%   predExprs(:,2:end) = model.comp{j}.ypred;
%   predFdelay = zeros(size(predF));
%   predFdelay(1) = predF(1);
%   predFdelay(2:end) = predF(1:end-1);
%   df = (predF - predFdelay)/0.1;
%   
%   predDi = (df + model.comp{j}.delta*predF)/ ...
%       model.comp{j}.sigma;
%   scaleExprs = sqrt(var(predExprs))./sqrt(var(y{j}));
%   meanExprs = mean(predExprs);
%   
%   model.comp{j}.sigma = model.comp{j}.sigma*scaleExprs(1);
%   model.comp{j}.B = (model.comp{j}.B-model.comp{j}.D.*meanExprs(2:end))./ ...
%       scaleExprs(2:end) + model.comp{j}.D.*mean(y{j}(2:end));
%   model.comp{j}.B = model.comp{j}.B./scaleExprs(2:end);
%   model.comp{j}.S = model.comp{j}.S./scaleExprs(2:end);
%   
%   model.comp{j}.kern.comp{2}.di_variance = model.comp{j}.sigma^2;
%   for i = 2:model.comp{j}.kern.numBlocks
%     model.comp{j}.kern.comp{i}.variance = model.comp{j}.S(i-1)^2;
%   end
%    
%   model.comp{j} = gpsimMapUpdateYpred(model.comp{j});
%   predExprs(:,2:end) = model.comp{j}.ypred;  
%   predExprs(:,1) = (df + model.comp{j}.delta*model.comp{j}.g)/ ...
%       model.comp{j}.sigma;  
   

  % Predicted Gene Expressions
  Kxx = multiKernCompute(model.comp{j}.kern, predt, model.comp{j}.t);
  meanPredX = reshape(ones(length(predt),1)*([0 model.comp{j}.B./ ...
                      model.comp{j}.D]), length(predt)*(numGenes+1), 1);
  predX = meanPredX + real(Kxx*model.comp{j}.invK*(model.comp{j}.y-ymean));
  varX = real(kernDiagCompute(model.comp{j}.kern, predt) - sum(Kxx'.* ...
                           (model.comp{j}.invK*Kxx'), 1)');

  % Take out predictions at data points.
  % Use them to get the scale for the other data.
  numData = length(predX)/(numGenes+1);
  predExprs = reshape(predX, numData, numGenes+1);
  meanExprs = ones(numData, 1)*mean(predExprs);
  scaleExprs = ones(numData, 1)*(sqrt(var(predExprs))./sqrt(var(y{j})));
  % Driving input can only be adjusted by the scale constant.
%  predExprs(:,1) = predExprs(:,1)./scaleExprs(:,1);
  % predictions of other genes can be generally normalised.
%   predExprs(:,2:end) =  ones(numData, 1)*mean(y{j}(:,2:end)) + ...
%       (predExprs(:,2:end) - meanExprs(:,2:end))./scaleExprs(:,2:end);
  predExprs(end-length(model.comp{j}.t)+1:end,:) = [];
  varExprs = reshape(varX, numData, numGenes+1);
  varExprs = varExprs./scaleExprs./scaleExprs;
  varExprs(end-length(model.comp{j}.t)+1:end,:) = []; 
  predF(end-length(model.comp{j}.t)+1:end) = [];
  varF(end-length(model.comp{j}.t)+1:end) = [];
  predt(end-length(model.comp{j}.t)+1:end) = [];
  
  yscale = ones(length(times), 1)*sqrt(var(y{j}));
  ymean = ones(length(times),1)*mean(y{j});
  ynormal = ymean + (y{j}-ymean)./yscale;
  yvarNormal = yvar{j}./yscale./yscale;

  scalePred = sqrt(var(predExprs));
  
  figure;
%  subplot(length(model.comp), 1, j);
  lin = plot(predt, predF, '-');
  hold on,
  bh = plot(predt, predF + 2*sqrt(varF), '--');
  bh = [bh plot(predt, predF - 2*sqrt(varF), '--')];
  hold off
%  ylabel(['Replica' num2str(j)]);
  title('Inferred Mef2 Protein', 'fontsize', 20);
  set(bh, 'lineWidth', 2);
  set(lin, 'lineWidth', 3);
  %set(lin, 'markersize', 20);
  set(gca, 'fontname', 'arial', 'fontsize', 24, 'xlim', [min(model.comp{j}.t) ...
                                                    max(model.comp{j} ...
                                                    .t)]);
  set(gca, 'ylim', [-0.1 0.4]);
  if saveFigures==1
    saveFileName = [fileName 'TF_profile_Rep' num2str(j)];
    print('-dpng', ['./results/' saveFileName]);
    pos = get(gcf, 'paperposition');
    origpos = pos;
    pos(3) = pos(3);
    pos(4) = pos(4);
    set(gcf, 'paperposition', pos);
    lineWidth = get(gca, 'lineWidth');
    set(gca, 'lineWidth', lineWidth*2);
    print('-deps', ['./results/' saveFileName]);  
    set(gca, 'lineWidth', lineWidth);
    set(gcf, 'paperposition', origpos);
  end
  
  for index = 1:numGenes+1
    figure;
    lin = plot(predt, predExprs(:,index), '-');
    hold on,
    bh = plot(predt, predExprs(:,index)+2*real(sqrt(varExprs(:,index))), '--');
    bh = [bh plot(predt, predExprs(:,index)-2*real(sqrt(varExprs(:,index))), '--')];
    lin = [lin plot(times, y{j}(:,index), 'rx')];
    lin1 = errorbar(times, y{j}(:,index), 2*sqrt(yvar{j}(:,index)), ...
                    'rx');
    if index == 1
      titleText = ['Driving Input mRNA'];      
%      lin = [lin plot(predt, predDi, 'm-')];     
    else
      titleText = ['Gene ' genenames{index} ' mRNA'];
    end
    title(titleText, 'fontsize', 20);
    set(bh, 'lineWidth', 3);
    set(lin, 'lineWidth', 4);
    set(lin, 'markersize', 20);
    set(lin1, 'lineWidth', 2);  
    set(gca, 'fontname', 'arial', 'fontsize', 24, 'xlim', [min(predt) ...
                        max(predt)]);
    
    if saveFigures
      saveFileName = [fileName '_ExprsProfile_Rep' num2str(j) '_Gene' num2str(index)];
      print('-deps', ['./results/' saveFileName]);
      pos = get(gcf, 'paperposition');
      origpos = pos;
      pos(3) = pos(3);
      pos(4) = pos(4);
      set(gcf, 'paperposition', pos);
      lineWidth = get(gca, 'lineWidth');
      set(gca, 'lineWidth', lineWidth*2);
      print('-dpng', ['./results/' saveFileName]);
      set(gca, 'lineWidth', lineWidth);
      set(gcf, 'paperposition', origpos);
    end    
  end
    
end

