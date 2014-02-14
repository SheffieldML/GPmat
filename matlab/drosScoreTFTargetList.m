function results = drosScoreTFTargetList(tf, listfile, savepath, modulus, reminder, subsample, noiseModel),

% DROSSCORETFTARGETLIST Score a list of candidate targets.
% FORMAT
% DESC Score a list of candidate targets of a TF.
% Suitable for non-interactive use, e.g. on a cluster when saving
% results to a file. In this mode, interrupted runs can be continued
% by rerunning with same parameters.
% ARG tf : TF symbol, should be in {'bap', 'bin', 'mef2', 'tin', 'twi'}
% ARG listfile : If a string: full path to a list of targets to test.
% If a cell array: a list of targets.
% ARG savepath : Full path of a directory where results are saved
% (set to empty to return a structure instead)
% ARG modulus : Number of parallel jobs in the set this one belongs, see below
% ARG reminder : Index of this parallel job in the set 0, ..., modulus-1
% This run will only look at genes where mod(index, modulus) == reminder
% Use modulus=1, reminder=0 to disable parallelisation.
% ARG subsample : Time point indices to use (default: [] = all)
% ARG noiseModel : Noise model to use (default: drosGpdisimLearn default)
% RETURN results : a structure with the fields:
% ll : a list of model log-likelihoods
% params : a cell array of model parameters
% targets : a cell array of target labels
% subsample : indices of data used (empty = all)
%
% Additionally if savepath is set, writes file
% 'dros_gpdisim_{tf}_list_{basename(listfile)}_results.mat' (non-parallel) or
% 'dros_gpdisim_{tf}_list_{basename(listfile)}_{modulus}_{reminder}_results.mat'
% under savepath. The file contains the above fields as separate variables.
%
% SEEALSO : drosLoadData, drosGpdisimLearn, drosScoreTFTargetListMT
%
% COPYRIGHT : Antti Honkela, 2009

% SHEFFIELDML

if nargin < 4,
  modulus = 1;
  reminder = 0;
end

if nargin < 6,
  subsample = [];
end

drosLoadData;

if ~any(strcmp(drosTF.names, tf)),
  error('Unknown TF');
end

if ~iscell(listfile),
  [foo1, infile, foo2, foo3] = fileparts(listfile);
else
  infile = 'cellinput';
end

if ~isempty(savepath),
  if modulus == 1,
    fname = sprintf('%s/dros_gpdisim_%s_list_%s_results.mat', savepath, tf, infile);
  else
    fname = sprintf('%s/dros_gpdisim_%s_list_%s_%d_%d_results.mat', savepath, tf, infile, modulus, reminder);
  end
  if exist(fname, 'file'),
    load(fname);
  else
    basefname = sprintf('%s/dros_gpdisim_%s_list_%s_results.mat', savepath, tf, infile);
    if exist(fname, 'file'),
      load(fname);
    else
      ll = [];
      params = {};
      if ~iscell(listfile),
	targets = importdata(listfile);
      else
	targets = listfile;
      end
    end
  end
else
  fname = [];
  ll = [];
  params = {};
  if ~iscell(listfile),
    targets = importdata(listfile);
  else
    targets = listfile;
  end
end

gpdisimArgs = struct('randomize', {0}, ...
		     'subsample', {subsample})

if nargin > 6,
  gpdisimArgs.noiseModel = noiseModel;
else
  noiseModel = [];
end

for l=1:length(targets),
  if ((l <= length(ll)) && isfinite(ll(l)) && ~isempty(params{l})) || (mod(l, modulus) ~= reminder),
    continue;
  end
  if any(strcmp(drosexp.probes, targets{l})),
    fprintf('Scoring %s...\n', targets{l});
    try,
      [ll(l), params{l}] = learn_it(drosexp, drosTF, tf, targets(l), gpdisimArgs);
    catch
      % Something went wrong, let's try again with a different initialisation
      for k=1:10,
	fprintf('Retrying gene %s\n', targets{l});
	failed = 0;
	try
	  gpdisimArgs.randomize = 1;
	  [ll(l), params{l}] = learn_it(drosexp, drosTF, tf, targets(l), gpdisimArgs);
	catch
	  failed = 1;
	end
	gpdisimArgs.randomize = 0;
	if ~failed,
	  break;
	end
      end
    end
  else
    ll(l) = -Inf;
    params{l} = [];
  end
  if ~isempty(savepath),
    save(fname, 'll', 'targets', 'params', 'subsample', 'noiseModel');
  end
end

results = struct('ll', {ll}, 'targets', {targets}, ...
		 'params', {params}, 'subsample', {subsample}, ...
		 'noiseModel', {noiseModel});



function [ll, params] = learn_it(drosexp, drosTF, tf, targets, args),

m = drosGpdisimLearn(drosexp, drosTF, tf, targets, args);
ll = modelLogLikelihood(m);
if ~isfinite(ll),
  error('non-finte log-likelihood');
end
params = modelExtractParam(m);
