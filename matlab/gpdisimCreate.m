function model = gpdisimCreate(numGenes, numProteins, times, geneVals, ...
			       geneVars, options, annotation)

% GPDISIMCREATE Create a GPDISIM model.
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
% ARG geneVars : the varuabces of each gene at the different time points.
% ARG options : options structure, the default options can be
% generated using gpsimOptions.
% ARG annotation : annotation for the data (gene names, etc.) that
% is stored with the model. (Optional)
% RETURN model : model structure containing default
% parameterisation.
%
% SEEALSO : modelCreate, gpsimOptions
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Antti Honkela, 2007

% GPSIM

if any(size(geneVars)~=size(geneVals))
  error('The gene variances have a different size matrix to the gene values.');
end

if(numGenes ~= (size(geneVals, 2) - 1))
  error('The number of genes given does not match the dimension of the gene values given.')
end

if(size(times, 1) ~= size(geneVals, 1))
  error('The number of time points given does not match the number of gene values given')
end

model.type = 'gpdisim';

kernType1{1} = 'multi';
kernType2{1} = 'multi';
tieWidth = [1]; % These are the indices of the inverse widths which
                % need to be constrained to be equal.
tieRBFVariance = [2];
kernType1{2} = 'rbf';
for i = 1:numGenes
  kernType1{i+2} = 'disim';
  if i==1
    tieDelta = [3];
    tieWidth = [tieWidth, 4];
    tieSigma = [5];
    tieRBFVariance = [tieRBFVariance, 8];
  end
  if i>1
    tieDelta = [tieDelta tieDelta(end)+6];
    tieWidth = [tieWidth tieWidth(end)+6];
    tieSigma = [tieSigma tieSigma(end)+6];
    tieRBFVariance = [tieRBFVariance tieRBFVariance(end)+6];
  end
end
tieParam = {tieDelta, tieWidth, tieSigma, tieRBFVariance};

model.y = geneVals(:);
model.yvar = geneVars(:);

model.includeNoise = options.includeNoise;

if model.includeNoise
  model.yvar = zeros(size(geneVars(:)));
else
  model.yvar = geneVars(:);
end

% Check if we have a noise term.
if model.includeNoise
  % Create a new multi kernel to contain the noise term.
  kernType2{1} = 'multi';

  % Set the new multi kernel to just contain 'white' kernels.
  for i = 1:numGenes+1
    kernType2{i+1} = 'white';
  end
  if isfield(options, 'singleNoise') & options.singleNoise
    tieParam{5} = (2+6*numGenes + 1):(2+6*numGenes + numGenes+1);
  end
  
  % Now create model with a 'cmpnd' (compound) kernel build from two
  % multi-kernels. The first multi-kernel is the sim-sim one the next
  % multi-kernel is the white-white one. 
  model.kern = kernCreate(times, {'cmpnd', kernType1, kernType2});
  simMultiKern = 'model.kern.comp{1}';
else
  model.kern = kernCreate(times, kernType1);
  simMultiKern = 'model.kern';
end

% This is if we need to place priors on parameters ...
if isfield(options, 'addPriors') && options.addPriors,
  for i = 1:length(model.kern.numBlocks)
    % Priors on the sim kernels.
    eval([simMultiKern '.comp{i}.priors = priorCreate(''gamma'');']);
    eval([simMultiKern '.comp{i}.priors.a = 1;']);
    eval([simMultiKern '.comp{i}.priors.b = 1;']);
    %model.kern.comp{i}.priors = priorCreate('gamma');
    %model.kern.comp{i}.priors.a = 1;
    %model.kern.comp{i}.priors.b = 1;
    if i == 1
      % For first kernel place prior on inverse width.
      % model.kern.comp{i}.priors.index = [1 2];
      eval([simMultiKern '.comp{i}.priors.index = [1 2];']);
    elseif i == 2
      %model.kern.comp{i}.priors.index = [1 3 4 5];
      eval([simMultiKern '.comp{i}.priors.index = [1 3 4 5];']);
    else
      % For other kernels don't place prior on inverse width --- as
      % they are all tied together and it will be counted multiple
      % times.
      %model.kern.comp{i}.priors.index = [4 5];
      eval([simMultiKern '.comp{i}.priors.index = [4 5];']);
    end
  end

  % Prior on the b values.
  model.bprior = priorCreate('gamma');
  model.bprior.a = 1;
  model.bprior.b = 1;
end

model.kern = modelTieParam(model.kern, tieParam);
if model.includeNoise,
  model.kern.comp{2}.comp{1}.variance = 1e-6;
end

% The decays and sensitivities are actually stored in the kernel.
% We'll put them here as well for convenience.
model.delta = 10;
model.sigma = 1;
for i = 2:model.kern.numBlocks
  eval([simMultiKern '.comp{i}.di_decay = model.delta;']);
  eval([simMultiKern '.comp{i}.di_variance = model.sigma^2;']);
  model.D(i-1) = eval([simMultiKern '.comp{i}.decay;']);
  model.S(i-1) = sqrt(eval([simMultiKern '.comp{i}.variance;']));
end

rand('seed',0);
model.numParams = numGenes + model.kern.nParams;
model.numGenes = numGenes;
model.mu = mean(geneVals(:, 2:end));
% model.B = model.D.*model.mu;
model.B = model.D.*geneVals(1, 2:end);
model.m = model.y;
model.t = times;

model.optimiser = options.optimiser;

if isfield(options, 'fix')
  model.fix = options.fix;
end

% The basal transcriptions rates must be postitive.
model.bTransform = optimiDefaultConstraint('positive');

if nargin > 6,
  model.annotation = annotation;
end

model.options = options;

% This forces kernel compute.
params = gpdisimExtractParam(model);
model = gpdisimExpandParam(model, params);

