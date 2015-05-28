% DEMBARENCOVARIATIONAL1 Optimise model using variational approximation with RBF kernel and exponential response.

% SHEFFIELDML

% OPTIMISATION OF PARAMETERS NOT YET IMPLEMENTED!!!

expNo = 1;
type = 'variational';

load demBarenco1;
origModel = model;
clear model

colordef white
[y, yvar, gene, times, scale, rawExp, rawVar] = gpsimLoadBarencoData;

numGenes = size(gene, 1);
% Get the default options structure.
options = gpsimMapOptions(numGenes);
options.kern = {'exp', 'rbf'};
options.meanFunction = 'expvarMean';
% Options for the expvarMean function are just the kernel.
options.meanFunctionOptions.kern = options.kern;
options.nonLinearity = 'linear';
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

for i = 1:3
  model.comp{i}.kern.argument.variance = 1.5;
  model.comp{i}.meanFunction.kern.argument.variance = 1.5;
end

paramvec{1} = gpsimMapExtractParam(model.comp{1}); %vector of gamma estimates

eta=0.02;
  
for rep=1:length(model.comp) %Work out likelihood gradient for each replicate    
  options = defaultOptions;
  options(1) = 1;
  model.comp{rep} = gpsimMapUpdateF(model.comp{rep}, options);
  ll(rep) = gpsimMapLogLikelihood(model.comp{rep});
  
  dg{rep} = gpsimMapLogLikeGradients(model.comp{rep});    
end 


fprintf('Log-likelihood %2.4f\t%2.4f\t%2.4f\n', ...
        ll(1), ll(2), ll(3));

type(1) = upper(type(1));  
save(['demBarenco' type num2str(expNo)], 'model', 'type', 'expNo')

dummyModel = model;
for rep = 1:length(dummyModel.comp)
  postMean = dummyModel.comp{rep}.f + modelOut(dummyModel.comp{rep}.meanFunction, ...
                                          dummyModel.comp{rep} ...
                                               .mapt);
  postMean(find(postMean<0)) = eps;
  postVar = dummyModel.comp{rep}.varf;
  
  % "Invert" log normal approximation.
  dummyModel.comp{rep}.varf = log(1+postVar./(postMean.*postMean));
  dummyModel.comp{rep}.f = log(postMean) - .5*dummyModel.comp{rep}.varf;
  dummyModel.comp{rep}.nonLinearity = 'exp';
end
gpsimMapBarencoResults(dummyModel, type, expNo)
