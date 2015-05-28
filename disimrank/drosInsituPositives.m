function v = drosInsituPositives(ids, drosinsitu, drosexp),

% DROSINSITUPOSITIVES Filter genes by looking for positive in-situ annotations.
% FORMAT
% DESC Filter a set of genes by retaining ones with positive in-situ
% annotations, keeping the order of the original list.
% ARG ids : Genes to look, expressed as indices in drosexp
% ARG drosinsitu : drosinsitu data set returned by drosLoadData
% ARG drosexp : drosexp data set returned by drosLoadData
% RETURN v : the filtered indices
%
% SEEALSO : drosLoadData
%
% COPYRIGHT : Antti Honkela, 2009

% SHEFFIELDML

isgenes = drosinsitu.genes(find(sum(drosinsitu.data, 2)));
I = drosFindGeneinds(drosexp, isgenes, 0);

[C, IA, IB] = intersect(ids, I);

v = ids(sort(IA));
