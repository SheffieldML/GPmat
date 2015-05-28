function model = drosGpdisimLearnWithFixedInit(drosexp, drosTF, tf, targets, initparams, fixcomps, randomize),

% DROSGPDISIMLEARNWITHFIXEDINIT Learns a GPDISIM model for given genes.
% FORMAT
% DESC Learns a GPDISIM model for given genes while fixing given parameters.
% ARG drosexp : the drosexp data structure from drosLoadData
% ARG drosTF : the drosTF data structure from drosLoadData
% ARG tf : TF symbol, should be in {'bap', 'bin', 'mef2', 'tin', 'twi'}
% ARG targets : a cell array of target gene probe names
% ARG initparams : initial parameters to fix (use NaN for free parameters)
% ARG fixcomps : indices of components whos parameters are known to be fixed
% (speeds up learning by not re-evaluating the kernel between these)
% ARG randomize : flag if random initialisation should be used 
% (optional, default = false)
% RETURN model : the learned GPDISIM model.
%
% SEEALSO : drosLoadData, drosScoreTFTargetListMT
%
% COPYRIGHT : Antti Honkela, 2009

% SHEFFIELDML

if nargin < 7,
  randomize = 0;
end

tfprobe = drosTF.probes.(tf);

if ~iscell(targets),
  targets = {targets};
end

genes = [tfprobe, targets(:)'];

genenames = genes;

[y, yvar, gene, times, scale, rawExp, rawVar] = gpdisimGetDrosData(drosexp, genes);

% Get the default options structure.
options = gpsimOptions;
options.includeNoise = 1;

options.fix(1).index = 4;
options.fix(1).value = expTransform(1, 'xtoa');

% Fix first output sensitivity to fix the scaling
I = find(~isnan(initparams));
for k=1:length(I),
  options.fix(k+1).index = I(k);
  options.fix(k+1).value = initparams(I(k));
end

options.fixBlocks = fixcomps;
% initialise the model.
model.type = 'cgpdisim';  % This new model type is a hack to run
                        % the model in a hierarchical manner.
                        % need to do this more elegantly later.
for i =1:3          %% 3 original
  model.comp{i} = gpdisimCreate(length(targets), 1, times, y{i}, yvar{i}, options, struct('tf', {tf}, 'genes', {genes}));
end

if randomize,
  a = modelExtractParam(model);
  a = randn(size(a));
  model = modelExpandParam(model, a);
end

% Learn the model.
model = modelOptimise(model, [], [], 1, 3000);
