function ids = drosMapGenesToIDs(genes),

% DROSMAPGENESTOIDS Map Drosophila FBgn gene identifiers to integers
% FORMAT
% DESC Map Drosophila FBgn gene identifiers to integers by dropping
% the FBgn prefix, i.e. "FBgn0012345" -> 12345
% ARG genes : A cell array of gene identifiers to map
% RETURN ids : An array of corresponding integer identifiers
%
% COPYRIGHT : Antti Honkela, 2009

% SHEFFIELDML

ids = zeros(size(genes));
for k=1:length(genes),
  if strcmp(genes{k}, 'NA'),
    ids(k) = NaN;
  else
    ids(k) = sscanf(genes{k}(5:end), '%d');
  end
end
