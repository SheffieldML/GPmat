function model = drosGpdisimLearn(drosexp, drosTF, tf, targets, randomize, addPriors, subsample),

% DROSGPDISIMLEARN Learns the single-target GPDISIM model for one gene.
% FORMAT
% DESC Learns the single-target GPDISIM model for one gene.
% ARG drosexp : the drosexp data structure from drosLoadData
% ARG drosTF : the drosTF data structure from drosLoadData
% ARG tf : TF symbol, should be in {'bap', 'bin', 'mef2', 'tin', 'twi'}
% ARG targets : a cell array of  target gene probe names
% ARG randomize : flag if random initialisation should be used 
% (optional, default = false)
% ARG addPriors : flag if priors should be added to parameters
% (optional, default = false)
% ARG subsample : which subset of time points to use
% (optional, default = [], meaning use all)
% RETURN model : the learned GPDISIM model.
%
% SEEALSO : drosLoadData, drosScoreTFTargetList
%
% COPYRIGHT : Antti Honkela, 2009

% DISIMRANK

if nargin < 5,
  randomize = 0;
end

if nargin < 6,
  addPriors = 0;
end

if nargin < 7,
  subsample = [];
end

tflabel = drosTF.labels(strcmp(tf, drosTF.names));

if ~iscell(targets),
  targets = {targets};
end

genes = [tflabel, targets(:)'];

genenames = genes;

[y, yvar, gene, times, scale, rawExp, rawVar] = gpdisimGetDrosData(drosexp, genes);

% Get the default options structure.
options = gpsimOptions;
options.includeNoise = 1;
% Fix one decay (from the fourth gene --- p21) to 0.8 hr^-1, and
% the corresponding sensitivity (see just after eqn 2 in the
% mathematical methods of Barenco et al.)
%options.fix(1).index = 2;
%options.fix(1).value = expTransform(1, 'xtoa');

% Fix sigma = 1 to normalise the scale of the protein
options.fix(1).index = 4;
options.fix(1).value = expTransform(1, 'xtoa');

if addPriors,
  options.addPriors = 1;
end

I = drosGetGeneinds(drosexp, genes);

% initialise the model.
model.type = 'cgpdisim';  % This new model type is a hack to run
                        % the model in a hierarchical manner.
                        % need to do this more elegantly later.
for i =1:3          %% 3 original
  if isempty(subsample),
    model.comp{i} = gpdisimCreate(length(targets), 1, times, y{i}, yvar{i}, options, struct('tf', {tf}, 'probes', {genes}, 'genes', {drosexp.fbgns(I)}, 'genenames', {drosexp.symbols(I)}));
  else
    model.comp{i} = gpdisimCreate(length(targets), 1, times(subsample), y{i}(subsample, :), yvar{i}(subsample, :), options, struct('tf', {tf}, 'probes', {genes}, 'genes', {drosexp.fbgns(I)}, 'genenames', {drosexp.symbols(I)}));
  end
end

if randomize,
  a = modelExtractParam(model);
  I = a==0;
  a = randn(size(a));
  a(I) = 0;
  model = modelExpandParam(model, a);
end

% Learn the model.
model = modelOptimise(model, [], [], 1, 3000);
