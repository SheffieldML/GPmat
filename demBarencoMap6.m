% DEMBARENCOMAP6 Optimise model using MAP approximation with MLP kernel and negLogLogit response.

% GPSIM

expNo = 6;
type = 'map';

load demBarenco1;
origModel = model;
clear model

colordef white
[y, yvar, gene, times, scale, rawExp, rawVar] = gpsimLoadBarencoData;

numGenes = size(gene, 1);
% Get the default options structure.
options = gpsimMapOptions(numGenes);
options.kern = 'mlp';
options.nonLinearity = 'negLogLogit';
options.includeNoise = 1;
options.intPoints = 161;

for i =1:3
  times = times;
  options.B = origModel.comp{1}.B;
  options.B = options.B;
  options.D = origModel.comp{1}.D;
  options.D = options.D;
  options.S = origModel.comp{1}.S;
  model.comp{i} = gpsimMapCreate(numGenes, 1, times, y{i}, yvar{i}, options);
  if strcmp(options.kern, 'mlp')
    model.comp{i}.kern.weightVariance = 30;
    model.comp{i}.kern.biasVariance = 1000;
    % This forces kernel recompute.
    params = gpsimMapExtractParam(model.comp{i});
    model.comp{i} = gpsimMapExpandParam(model.comp{i}, params);
  end

end

paramvec{1} = gpsimMapExtractParam(model.comp{1}); %vector of gamma estimates

eta=0.02;
iters = 50; %Number of optimisation iterations
for ii=1:iters %Start optimisation
  
  for rep=1:length(model.comp) %Work out likelihood gradient for each replicate    
    options = defaultOptions;
    options(1) = 1;
    model.comp{rep} = gpsimMapUpdateF(model.comp{rep}, options);
    ll(ii, rep) = gpsimMapLogLikelihood(model.comp{rep});
    
    dg{rep} = gpsimMapLogLikeGradients(model.comp{rep});    
  end 
  fprintf('Iteration %d, log-likelihood %2.4f\t%2.4f\t%2.4f\n', ...
          ii, ll(ii, 1), ll(ii, 2), ll(ii, 3));
  
  % Update kernel parameters by simple gradient ascent
  param = gpsimMapExtractParam(model.comp{rep});
  for i = 1:length(dg)
    param(1:end-1) = param(1:end-1) + eta*sum(dg{i}(1:end-1));
  end
  for rep = 1:length(model.comp)
    model.comp{rep} = gpsimMapExpandParam(model.comp{rep}, ...
                                          param);
  end
  paramvec{end+1}=param;
end 
type(1) = upper(type(1));  
save(['demBarenco' type num2str(expNo)], 'model', 'type', 'expNo')
gpsimMapBarencoResults(model, type, expNo)