function model = drosGpdisimLearn(drosexp, drosTF, tf, targets, varargin),

% DROSGPDISIMLEARN Learns a GPDISIM model for given genes.
% FORMAT
% DESC Learns a GPDISIM model for given genes.
% ARG drosexp : the drosexp data structure from drosLoadData
% ARG drosTF : the drosTF data structure from drosLoadData
% ARG tf : TF symbol, should be in {'bap', 'bin', 'mef2', 'tin', 'twi'}
% ARG targets : a cell array of target gene probe names
% ARG ... : Keyword, value pairs; following keywords are supported:
% 'randomize' : flag if random initialisation should be used (default = false)
% 'addPriors' : flag if priors should be added to parameters (default = false)
% 'subsample' : which subset of time points to use (default = [], meaning all)
% 'noiseModel' : noise model to use, possible values are
%   'probe' : use noise extracted from expression data by mmgMOS (default)
%   'shared' : use one global noise parameter
%   'individual' : use one noise parameter for each gene
% RETURN model : the learned GPDISIM model.
%
% SEEALSO : drosLoadData, drosScoreTFTargetList, drosScoreTFTargetListMT
%
% COPYRIGHT : Antti Honkela, 2009

% SHEFFIELDML

% Read the arguments
if (length(varargin) == 1) & (isstruct(varargin{1})),
  args = varargin{1};
elseif (mod(length(varargin), 2) ~= 0),
  if ~(isstruct(varargin{1})),
    error('Keyword arguments should appear in pairs');
  else
    args = varargin{1};
    for k=2:2:length(varargin),
      if ~ischar(varargin{k})
	error('Keyword argument names must be strings.');
      end
      args.(varargin{k}) = varargin{k+1};
    end
  end
else
  args = struct(varargin{:});
end

tfprobe = drosTF.probes.(tf);

if ~iscell(targets),
  targets = {targets};
end

genes = [tfprobe, targets(:)'];

genenames = genes;

% Get the default options structure.
options = gpsimOptions;
% Fix one decay (from the fourth gene --- p21) to 0.8 hr^-1, and
% the corresponding sensitivity (see just after eqn 2 in the
% mathematical methods of Barenco et al.)
%options.fix(1).index = 2;
%options.fix(1).value = expTransform(1, 'xtoa');

% Fix sigma = 1 to normalise the scale of the protein
options.fix(1).index = 4;
options.fix(1).value = expTransform(1, 'xtoa');

if isfield(args, 'addPriors') && args.addPriors,
  options.addPriors = 1;
end

if isfield(args, 'noiseModel'),
  switch args.noiseModel
   case 'probe',
    options.includeNoise = 0;
   case 'shared',
    options.includeNoise = 1;
    options.singleNoise = true;
   case 'individual',
    options.includeNoise = 1;
    options.singleNoise = false;
  end
else
  options.includeNoise = 0;
end

I = drosFindGeneinds(drosexp, genes, 0, 1);
if options.includeNoise,
  [y, yvar, gene, times, scale, rawExp, rawVar] = drosGetData(drosexp, I, 1);
else
  [y, yvar, gene, times, scale, rawExp, rawVar] = drosGetData(drosexp, I);
end

% initialise the model.
model.type = 'cgpdisim';  % This new model type is a hack to run
                        % the model in a hierarchical manner.
                        % need to do this more elegantly later.
for i =1:3          %% 3 original
  if isfield(args, 'subsample') && ~isempty(args.subsample),
    model.comp{i} = gpdisimCreate(length(targets), 1, times(args.subsample), y{i}(args.subsample, :), yvar{i}(args.subsample, :), options, struct('tf', {tf}, 'probes', {drosexp.probes(I)}, 'genes', {drosexp.genes(I)}, 'genenames', {drosexp.symbols(I)}));
  else
    model.comp{i} = gpdisimCreate(length(targets), 1, times, y{i}, yvar{i}, options, struct('tf', {tf}, 'probes', {drosexp.probes(I)}, 'genes', {drosexp.genes(I)}, 'genenames', {drosexp.symbols(I)}));
  end
end

if isfield(args, 'randomize') && args.randomize,
  a = modelExtractParam(model);
  I = a==0;
  a = randn(size(a));
  a(I) = 0;
  model = modelExpandParam(model, a);
end

% Learn the model.
model = modelOptimise(model, [], [], 1, 3000);
