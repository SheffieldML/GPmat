function r = drosCorrelationRank(drosexp, drosTF, tf),

% DROSCORRELATIONRANK Rank targets by expression profile correlation.
% FORMAT
% DESC Rank targets by expression profile correlation.
% ARG drosexp : the drosexp data structure from drosLoadData
% ARG drosTF : the drosTF data structure from drosLoadData
% ARG tf : TF symbol, should be in {'bap', 'bin', 'mef2', 'tin', 'twi'}
% RETURN ranking : indeces of genes in drosexp in ranked order
%
% SEEALSO : drosLoadData, drosScoreTFTargetList
%
% COPYRIGHT : Antti Honkela, 2009

% SHEFFIELDML

N = length(drosexp.probes);

medprofiles = mean(reshape(drosexp.pctiles(:, :, 3), [N, 12, 3]), 3);
I = strcmp(drosTF.probes.(tf), drosexp.probes);

sigma = zeros(1, N);
for k=1:N,
  v = corrcoef(medprofiles(I, :)', medprofiles(k, :)');
  sigma(k) = v(1, 2);
end

[foo, J] = sort(sigma, 'descend');
r = J(2:end);
