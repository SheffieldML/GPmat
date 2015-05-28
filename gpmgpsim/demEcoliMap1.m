% DEMECOLIMAP1 Optimise model using MAP approximation with MLP kernel and multiple repression response, using SCG optimisation and Ecoli data set.

% SHEFFIELDML

clear;
expNo = 4;
type = 'mapFullEcoliOptimInit';
rand('seed',2);

load demBarenco1;
origModel = model;
clear model;
clear y;
clear yvar;
saveFigures = false;

colordef white
[y, gene, times, scale, rawExp] = gpsimLoadEcoliFullData;
Nrep = length(y);
for i = 1:Nrep
  yvar{i} = zeros(size(y{i}));
end

numGenes = length(gene);
% Get the default options structure.
options = gpsimMapOptions(numGenes);
options.kern = {'translate', 'mlp'};    % 'mlp';
options.nonLinearity = 'repression';
options.includeNoise = 0;
options.startPoint = 0;                 % TF starts from 0;
options.endPoint = 60;
options.priorProteinTimes = []';
options.priorProtein = []';
options.includeNoise = 1;
options.bTransform = 'exp';
options.alphaTransform = [];
% Modify according to the range of the TF
options.intPoints = (options.endPoint-options.startPoint)/0.5 + 1 ;
options.gParam = ones(1,numGenes);
options.times = times;

for i =1:length(y)
    times = times;
    options = gpsimMapInitParam(y{i}, yvar{i}, options);
    
    if strcmp(options.kern{2}, 'mlp')
%       options.fix(1).index = 20;         % S of gene4
%       options.fix(1).value = expTransform(1.75, 'xtoa');
%       options.fix(2).index = 21;          % D of gene4
%       options.fix(2).value = expTransform(0.23, 'xtoa');
      options.fix(1).index = 3;  
      options.fix(1).value =  expTransform(1, 'xtoa');
%       options.fix(4).index = 22;  % alpha of gene 4
%       options.fix(4).value =  expTransform(0.01, 'xtoa');
    elseif strcmp(options.kern, 'rbf')
      options.fix(1).index = 19;         % S of gene4
      options.fix(1).value = expTransform(1.75, 'xtoa');
      options.fix(2).index = 20;          % D of gene4
      options.fix(2).value = expTransform(0.23, 'xtoa');
      options.fix(3).index = 2;  
      options.fix(3).value =  expTransform(1, 'xtoa');
      options.fix(4).index = 21;  % alpha of gene 4
      options.fix(4).value =  expTransform(0.01, 'xtoa');
    end        
    
    model.comp{i} = gpsimMapCreate(numGenes, 1, times, y{i}, yvar{i}, options);
    if strcmp(options.kern{2}, 'mlp')
        model.comp{i}.kern.comp{1}.weightVariance = 30;

        model.comp{i}.kern.comp{1}.biasVariance = 20;
        model.comp{i}.kern.centre = -15;
        % This forces kernel recompute.
        params = gpsimMapExtractParam(model.comp{i});
        model.comp{i} = gpsimMapExpandParam(model.comp{i}, params);
    elseif strcmp(options.kern, 'rbf')
      model.comp{i}.kern.inverseWidth = 0.01;
      params = gpsimMapExtractParam(model.comp{i});
      model.comp{i} = gpsimMapExpandParam(model.comp{i}, params);
    end
end

param = 0;

for i = 1:Nrep
  

    option = defaultOptions;
    model.comp{i} = gpsimMapUpdateF(model.comp{i}, option);
   
    % modify the initialisations
    fmean =  mean(model.comp{i}.f);
%/~
    %     fsd = sqrt(var(model.comp{i}.f));
%    model.comp{i}.f = model.comp{i}.f - fmean;
%     model.comp{i}.f = model.comp{i}.f/fsd;
%     model.comp{i}.S = model.comp{i}.S/exp(fmean);
%     model.comp{i}.gParam = model.comp{i}.gParam/exp(fmean);
%    ypredMean = mean(model.comp{i}.ypred);
%~/
    ypredScale = sqrt(var(model.comp{i}.ypred(model.comp{i}.times_index,:)));
    model.comp{i}.B = (model.comp{i}.B - model.comp{i}.D.* ...
                       (ypredMean - ypredScale.*mean(model.comp{i}.y)))./ypredScale;
    
    if isfield(model.comp{i}, 'alpha')
      model.comp{i}.alpha = model.comp{i}.alpha./ypredScale;
    end
    for j = 1:numGenes
      if model.comp{i}.B(j) < 0
        model.comp{i}.alpha(j) =  model.comp{i}.alpha(j) + model.comp{i}.B(j);
        model.comp{i}.B(j) = 1e-6;
      end
    end
    model.comp{i}.S = model.comp{i}.S./ypredScale; 
    f = gpsimMapFunctionalExtractParam(model.comp{i});
    model.comp{i} = gpsimMapFunctionalExpandParam(model.comp{i}, f);
    
    [paramvec{i} paramNames] = gpsimMapExtractParam(model.comp{i}); %vector of gamma estimates
    param = param + paramvec{i};    
    
end

param = param/Nrep;
ll0 = gpsimMapObjective(param, model)
iters = 300;

option = optOptions;
option(14) = iters;
option(9) = 0;
option(1) = 1;

newparam = scg('gpsimMapObjective', param,  option, ...
               'gpsimMapGradients', model);
%/~ 
% [newparam, ll, index] = minimize(param', 'gpsimMapGradFuncWrapper', iters, mo% del);
%~/
for rep = 1:Nrep
    option = defaultOptions; 
    model.comp{rep} = gpsimMapExpandParam(model.comp{rep}, newparam);   
    model.comp{rep} = gpsimMapUpdateF(model.comp{rep}, option);
    model.comp{rep}.yvar = ones(length(times), 1)* ...
        model.comp{rep}.noiseVar;
    model.comp{rep} = gpsimMapUpdateYpredVar(model.comp{rep});
end

type(1) = upper(type(1));  
save(['dem' type num2str(expNo)], 'model', 'type', 'expNo', 'scale')
gpsimMapEcoliResults(model, type, expNo, scale, saveFigures)
