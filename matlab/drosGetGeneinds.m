function I = drosGetGeneinds(drosdata, genes),
% DROSGETGENEINDS returns indices of given genes in the given data
%
% Usage:
%   I = drosGetGeneinds(drosdata, genes)
% COPYRIGHT : Antti Honkela, 2007
  
% SHEFFIELDML

I = [];

for k=1:length(genes),
  I = [I, find(strcmp(genes{k}, drosdata.genes))];
end
