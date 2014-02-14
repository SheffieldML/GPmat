function r = drosMergeResults(rr),

% DROSMERGERESULTS Merge split results
% FORMAT
% DESC Merge split results
% ARG rr : cell array of result structures
% RETURN r : one merged result structure
%
% SEEALSO : drosScoreTFTargetList, drosScoreTFTargetListMT
%
% COPYRIGHT : Antti Honkela, 2009

% SHEFFIELDML

N = length(rr);

r = rr{1};
I = ~(r.ll == 0);
r.ll(isnan(r.ll)) = -Inf;

for k=2:N,
  rr{k}.ll(isnan(rr{k}.ll)) = -Inf;
  J = ~(rr{k}.ll == 0);
  L = max(length(I), length(J));
  I(end+1:L) = false;
  J(end+1:L) = false;
  if any(r.ll(I & J) ~= rr{k}.ll(I & J)) || ~all(strcmp(r.targets, rr{k}.targets)),
    error('drosMergeResults: Components are inconsistent');
  end
  r.ll(J) = rr{k}.ll(J);
  rr{k}.params(length(rr{k}.params)+1:length(J)) = {[]};
  r.params(J) = rr{k}.params(J);
  I = ~(r.ll == 0);
end
