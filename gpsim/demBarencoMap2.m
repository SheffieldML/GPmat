% DEMBARENCOMAP2 Optimise model using MAP approximation with MLP kernel and exp response, using SCG optimisation and new PUMA processed data.

% SHEFFIELDML
clear,clc
expNo = 2;
type = 'mapMLPExp';
rand('seed',1);

load demBarenco1;
origModel = model;
clear model
saveFigures = true;

colordef white
[y, yvar, gene, times, scale, rawExp, rawVar] = gpsimLoadBarencoPUMAData;

numGenes = size(gene, 1);
% Get the default options structure.
options = gpsimMapOptions(numGenes);
options.kern = 'mlp';
options.nonLinearity = 'exp';
options.includeNoise = 0;
options.startPoint = 0;                 % TF starts from 0;
options.endPoint = 12;
options.priorProteinTimes = [0]';
options.priorProtein = [-10]';
% Modify according to the range of the TF
options.intPoints = (options.endPoint-options.startPoint)/0.1 + 1 ;
options.gParam = [];
options.ngParam = 0;

options.fix(1).index = 14;         % S of p21
options.fix(1).value = expTransform(1, 'xtoa');
options.fix(2).index = 15;          % D of p21
options.fix(2).value = expTransform(0.8, 'xtoa');
options.fix(3).index = 3;
options.fix(3).value =  expTransform(1, 'xtoa');
% options.fix(2).index = 1;
% options.fix(2).value =  expTransform(1, 'xtoa');;

for i =1:3
    times = times;
%    options.B = origModel.comp{1}.B;
%    options.D = origModel.comp{1}.D; % [0.6417 0.4094 0.4894 0.8000 0.4992];
%    options.S = origModel.comp{1}.S; % [0.4206 0.3255 0.1486 1 0.1870];
    options.S = ones(1, 5);    
    options.D = rand(1, 5);
    mu = mean(y{i}, 1);
    options.B = options.D.*mu;
    model.comp{i} = gpsimMapCreate(numGenes, 1, times, y{i}, yvar{i}, options);
    if strcmp(options.kern, 'mlp')
        model.comp{i}.kern.weightVariance = 30;
        model.comp{i}.kern.biasVariance = 1000;
        % This forces kernel recompute.
        params = gpsimMapExtractParam(model.comp{i});
        model.comp{i} = gpsimMapExpandParam(model.comp{i}, params);
    end

end

param = 0;
Nrep = length(model.comp);

for i = 1:Nrep
    paramvec{i} = gpsimMapExtractParam(model.comp{i}); %vector of gamma estimates
    param = param + paramvec{i};
end

param = param/Nrep;

iters = 300;

options = optOptions;
options(14) = iters;
options(9) = 0;
options(1) = 1;

newparam = scg('gpsimMapObjective', param,  options, ...
               'gpsimMapGradients', model);

% [newparam, ll, index] = minimize(param', 'gpsimMapGradFuncWrapper', iters, model);

for rep = 1:Nrep
    options = defaultOptions; 
    model.comp{rep} = gpsimMapExpandParam(model.comp{rep}, newparam);   
    model.comp{rep} = gpsimMapUpdateF(model.comp{rep}, options);
    model.comp{rep} = gpsimMapUpdateYpredVar(model.comp{rep});
end

type(1) = upper(type(1));  
save(['demBarenco' type num2str(expNo)], 'model', 'type', 'expNo', 'scale')
gpsimMapBarencoResults(model, type, expNo, saveFigures, scale)
