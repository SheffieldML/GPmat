function r = drosComputeChipDistanceCurve(chipdata, drosexp, drosinsitu, rankings, t, tf, filter);

% DROSCOMPUTECHIPDISTANCECURVE Compute the curve of distances to ChIP binding sites
% FORMAT
% DESC Compute the curve of distances to ChIP binding sites
% ARG chipdistances : the chipdistances data structure from drosLoadData
% ARG drosexp : the drosexp data structure from drosLoadData
% ARG drosinsitu : the drosinsitu data structure from drosLoadData
% ARG rankings : a cell array of R rankings, each of which is an array
% of indices of genes in drosexp in order of descending preference
% ARG t : the n in top-n to consider
% ARG tf : TF symbol, should be in {'bap', 'bin', 'mef2', 'tin', 'twi'}
% ARG filter : a flag whether to filter the results by positive in-situs
% (default: []=false)
% RETURN r : An Nx(R+1) array of frequencies of genes of binding within
% certain range, where the range (first index) ranges as 10.^[1:.1:6].
% The second index is over rankings, with the (R+1)'st one being random.
%
% SEEALSO : drosLoadData, demRunRankings, drosPlotChipDistances
%
% COPYRIGHT : Antti Honkela, 2009

% SHEFFIELDML

if nargin < 5,
  t = 100;
end

if nargin < 6,
  tf = 'twi';
end

if nargin < 7,
  filter = [];
end

chip_I = drosFindGeneinds(chipdata, drosexp.geneids, 1);
chip_J = strmatch(tf, chipdata.labels);

thresholds = 10.^[1:.1:6];
r = zeros(length(thresholds), length(rankings)+1);

if ~isempty(filter),
  for k=1:length(rankings),
    rankings{k} = drosInsituPositives(rankings{k}, drosinsitu, drosexp);
  end
end

validation = NaN * ones(size(drosexp.genes));
validation(chip_I~=0) = 0;

if ~isempty(filter),
  M = length(drosRemoveDuplicateGenes(drosexp, drosInsituPositives(find(~isnan(validation)), drosinsitu, drosexp)));
  M0 = length(drosRemoveDuplicateGenes(drosexp, find(~isnan(validation))));
else
  M = length(drosRemoveDuplicateGenes(drosexp, find(~isnan(validation))));
end

for m=1:length(thresholds),
  fprintf('Iteration %d/%d...\n', m, length(thresholds));
  threshold = thresholds(m);
  validation(chip_I~=0) = chipdata.data(chip_I(chip_I~=0), chip_J) < threshold;

  val2 = validation;
  val2(isnan(val2)) = 0;
  if ~isempty(filter),
    K = length(drosRemoveDuplicateGenes(drosexp, drosInsituPositives(find(val2), drosinsitu, drosexp)));
    K0 = length(drosRemoveDuplicateGenes(drosexp, find(val2)));
  else
    K = length(drosRemoveDuplicateGenes(drosexp, find(val2)));
  end

  for k=1:length(rankings),
    myranking = rankings{k};
    for l=1:length(t),
      if t(l) <= length(myranking),
	val = validation(myranking(1:t(l)));
	r(m, k) = nanmean(val);
	%pvals(l, k) = 1 - hygecdf(nansum(val)-1, M, K, sum(~isnan(val)));
      else
	r(m, k) = 0;
	%pvals(l, k) = 1;
      end
    end
  end
  r(m, k+1) = K/M;
end
