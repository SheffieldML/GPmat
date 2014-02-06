function model = gpsimCreate(numGenes, numProteins, times, geneVals, ...
                             geneVars, options, annotation)

% GPSIMCREATE Create a GPSIM model.
% The GPSIM model is a model for estimating the protein
% concentration in a small gene network where several genes are
% governed by one protein. The model is based on Gaussian processes
% and simple linear differential equations of the form
%
% dx(t)/dt = B + Cf(t) - Dx(t)
%
% where x(t) is a given genes concentration and f(t) is the protein
% concentration. 
%
% FORMAT
% DESC creates a model for single input motifs with Gaussian
% processes.
% ARG numGenes : number of genes to be modelled in the system.
% ARG numProteins : number of proteins to be modelled in the
% system.
% ARG times : the time points where the data is to be modelled.
% ARG geneVals : the values of each gene at the different time points.
% ARG geneVars : the variances of each gene at the different time points.
% ARG options : options structure, the default options can be
% generated using gpsimOptions.
% ARG annotation : annotation for the data (gene names, etc.) that
% is stored with the model. (Optional)
% RETURN model : model structure containing default
% parameterisation.
%
% SEEALSO : modelCreate, gpsimOptions
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2007
%
% MODIFIED : Pei Gao, 2008
% MODIFIED : Antti Honkela, 2008, 2009

% SHEFFIELDML

if any(size(geneVars)~=size(geneVals))
  error('The gene variances have a different size matrix to the gene values.');
end

if(numGenes ~= size(geneVals, 2))
  error('The number of genes given does not match the dimension of the gene values given.')
end

if(size(times, 1) ~= size(geneVals, 1))
  error('The number of time points given does not match the number of gene values given')
end

model.type = 'gpsim';

kernType1{1} = 'multi';

if isfield(options, 'proteinPrior') && ~isempty(options.proteinPrior)
  model.proteinPrior = options.proteinPrior;
  kernType1{2} = 'rbf';
  if isfield(options, 'proteinPriorTimes')
    timesCell{1} = options.proteinPriorTimes;
  else
    timesCell{1} = times;
  end
  for i = 1:numGenes
    kernType1{i+2} = 'sim';
    timesCell{i+1} = times; 
  end  
  model.timesCell = timesCell;
else
  timesCell = times;                     % Non-cell structure in this case
  for i = 1:numGenes
    kernType1{i+1} = 'sim';
  end
end
tieParam = {'inverse width'};

if isfield(options, 'fixBlocks') && ~isempty(options.fixBlocks),
  kernType1 = {'parametric', struct('fixBlocks', {options.fixBlocks}), kernType1};
end

model.y = geneVals(:);

model.includeNoise = options.includeNoise;

% if model.includeNoise
%   model.yvar = zeros(size(geneVars(:)));
% else
  model.yvar = geneVars(:);
% end

% Check if we have a noise term.
if model.includeNoise
  % Create a new multi kernel to contain the noise term.
  kernType2{1} = 'multi';

  % Set the new multi kernel to just contain 'white' kernels.
  if isfield(model, 'proteinPrior') && ~isempty(model.proteinPrior)
    kernType2{2}='whitefixed';
    for i = 2:(numGenes+1)
      kernType2{i+1} = 'white';
    end    
  else
    for i = 1:numGenes
      kernType2{i+1} = 'white';
    end
  end
  if isfield(options, 'singleNoise') & options.singleNoise
    tieParam{2} = 'white . variance';
  end
            
  
  % Now create model with a 'cmpnd' (compound) kernel build from two
  % multi-kernels. The first multi-kernel is the sim-sim one the next
  % multi-kernel is the white-white one. 
  model.kern = kernCreate(timesCell, {'cmpnd', kernType1, kernType2});
else
  model.kern = kernCreate(timesCell, kernType1);
end

% This is if we need to place priors on parameters ...
if isfield(options, 'addPriors') && options.addPriors,
  for i = 1:length(model.kern.numBlocks)
    % Priors on the sim kernels.
    model.kern.comp{i}.priors = priorCreate('gamma');
    model.kern.comp{i}.priors.a = 1;
    model.kern.comp{i}.priors.b = 1;
    if i == 1
      % For first kernel place prior on inverse width.
      model.kern.comp{i}.priors.index = [1 2 3];
    else
      % For other kernels don't place prior on inverse width --- as
      % they are all tied together and it will be counted multiple
      % times.
      model.kern.comp{i}.priors.index = [1 3];
    end
  end

  % Prior on the b values.
  model.bprior = priorCreate('gamma');
  model.bprior.a = 1;
  model.bprior.b = 1;
end

model.kern = modelTieParam(model.kern, tieParam);
model.kern.comp{2}.comp{1}.variance = 1e-6;

% The decays and sensitivities are actually stored in the kernel.
% We'll put them here as well for convenience.
if isfield(model, 'proteinPrior') && ~isempty(model.proteinPrior)
  for i = 2:model.kern.numBlocks
    if model.includeNoise
      model.D(i-1) = model.kern.comp{1}.comp{i}.decay;
      model.S(i-1) = sqrt(model.kern.comp{1}.comp{i}.variance);
    else
      model.D(i-1) = model.kern.comp{i}.decay;
      model.S(i-1) = sqrt(model.kern.comp{i}.variance);
    end
  end  
else
  for i = 1:model.kern.numBlocks
    if model.includeNoise
      model.D(i) = model.kern.comp{1}.comp{i}.decay;
      model.S(i) = sqrt(model.kern.comp{1}.comp{i}.variance);
    else
      model.D(i) = model.kern.comp{i}.decay;
      model.S(i) = sqrt(model.kern.comp{i}.variance);
    end
  end 
end

model.numParams = numGenes + model.kern.nParams;
model.numGenes = numGenes;
model.mu = mean(geneVals);
model.B = model.D.*model.mu;

if isfield(model, 'proteinPrior') && ~isempty(model.proteinPrior)
  dim = size(model.proteinPrior, 1) + size(model.y, 1);
  model.m = [model.proteinPrior; model.y];
else
  model.m = model.y;
  model.t = times;
end

model.optimiser = options.optimiser;

if isfield(options, 'fix')
  model.fix = options.fix;
end

% The basal transcriptions rates must be postitive.
model.bTransform = optimiDefaultConstraint('positive');

if nargin > 6,
  model.annotation = annotation;
end

% This forces kernel compute.
params = gpsimExtractParam(model);
model = gpsimExpandParam(model, params);

