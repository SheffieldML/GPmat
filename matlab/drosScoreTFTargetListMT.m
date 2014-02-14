function results = drosScoreTFTargetListMT(tf, listfile, basefile, savepath, modulus, reminder),

% DROSSCORETFTARGETLISTMT Score a list of candidate targets using multiple-target models.
% FORMAT
% DESC Score a list of candidate targets of a TF using multiple-target
% models and write the results to a file. Meant for non-interactive
% use, e.g. on a cluster. Interrupted runs can be continued by rerunning
% with the same parameters.
% ARG tf : TF symbol, should be in {'bap', 'bin', 'mef2', 'tin', 'twi'}
% ARG listfile : If a string: full path to a list of targets to test.
% If a cell array: a list of targets.
% ARG basefile : If a string: full path to a list of training set genes.
% If a cell array: a list of training set genes.
% ARG savepath : Full path of a directory where results are saved
% ARG modulus : Number of parallel jobs in the set this one belongs, see below
% ARG reminder : Index of this parallel job in the set 0, ..., modulus-1
% This run will only look at genes where mod(index, modulus) == reminder
% Use modulus=1, reminder=0 to disable parallelisation.
% RETURN results : a structure with the fields:
% ll : a list of model log-likelihoods
% params : a cell array of model parameters
% targets : a cell array of target labels
% basetargets : a cell array of training set target labels
%
% Additionally if savepath is set, writes file
% 'dros_gpdisimmt_{tf}_list_{basename(listfile)}_{basename(basefile)}_results.mat' (non-parallel) or
% 'dros_gpdisimmt_{tf}_list_{basename(listfile)}_{basename(basefile)}_{modulus}_{reminder}_results.mat'
% under savepath. The file contains the above fields as separate variables.
%
% SEEALSO : drosLoadData, drosGpdisimLearn, drosScoreTFTargetList
%
% COPYRIGHT : Antti Honkela, 2009

% SHEFFIELDML

if nargin < 5,
  modulus = 1;
  reminder = 0;
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
if ~iscell(basefile),
  [foo1, baseinfile, foo2, foo3] = fileparts(basefile);
else
  baseinfile = 'cellinput';
end

if ~isempty(savepath),
  if modulus==1,
    fname = sprintf('%s/dros_gpdisimmt_%s_list_%s_%s_results.mat', savepath, tf, infile, baseinfile);
  else
    fname = sprintf('%s/dros_gpdisimmt_%s_list_%s_%s_%d_%d_results.mat', savepath, tf, infile, baseinfile, modulus, reminder);
  end
  if exist(fname, 'file'),
    load(fname);
  else
    fname = sprintf('%s/dros_gpdisimmt_%s_list_%s_%s_results.mat', savepath, tf, infile, baseinfile);
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
  
  if modulus==1,
    fname = sprintf('%s/dros_gpdisimmt_%s_list_%s_%s_results.mat', savepath, tf, infile, baseinfile);
  else
    fname = sprintf('%s/dros_gpdisimmt_%s_list_%s_%s_%d_%d_results.mat', savepath, tf, infile, baseinfile, modulus, reminder);
  end
else
  if ~iscell(listfile),
    targets = importdata(listfile);
  else
    targets = listfile;
  end
end

if ~iscell(basefile),
  basetargets = importdata(basefile);
else
  basetargets = basefile;
end
basenum = length(basetargets);

m_base = drosGpdisimLearn(drosexp, drosTF, tf, basetargets);
a = modelExtractParam(m_base);

baseparams = NaN*ones(1, length(a) + 3);
baseparams(1:(2*basenum+4)) = a(1:(2*basenum+4));
t = 2*basenum + 5;
baseparams((t+2):(t+2+basenum-1)) = a(t:(t+basenum-1));

for l=1:length(targets),
  if ((l <= length(ll)) && isfinite(ll(l)) && ~isempty(params{l})) || (mod(l, modulus) ~= reminder),
    continue;
  end
  if any(strcmp(drosexp.probes, targets{l})),
    try,
      [ll(l), params{l}] = learn_it(drosexp, drosTF, tf, [basetargets; targets{l}], baseparams, 1:(basenum+1), 0);
    catch
      % Something went wrong, let's try again with a different initialisation
      for k=1:10,
	fprintf('Retrying gene %s\n', targets{l});
	failed = 0;
	try
	  [ll(l), params{l}] = learn_it(drosexp, drosTF, tf, [basetargets; targets{l}], baseparams, 1:(basenum+1), 1);
	catch
	  failed = 1;
	end
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
    save(fname, 'll', 'targets', 'params', 'basetargets');
  end
end

results = struct('ll', {ll}, 'targets', {targets}, ...
		 'params', {params}, 'basetargets', {basetargets});



function [ll, params] = learn_it(drosexp, drosTF, tf, targets, baseparams, fixcomps, randomize),

m = drosGpdisimLearnWithFixedInit(drosexp, drosTF, tf, targets, baseparams, fixcomps, randomize);
ll = modelLogLikelihood(m);
if ~isfinite(ll),
  error('non-finite log-likelihood');
end
params = modelExtractParam(m);
