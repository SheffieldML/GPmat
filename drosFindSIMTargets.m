function genes = drosFindSIMTargets(droschip, drosTF, tf),

% DROSFINDSIMTARGETS returns a list of targets only controlled by given tf
% (according to the given ChIP-chip data)
%
% Usage:
%   genes = drosFindSIMTargets(droschip, drosTF, tf)
% where tf is one of 'bap', 'bin', 'mef2', 'tin', 'twi'
%
% COPYRIGHT : Antti Honkela, 2007
  
% SHEFFIELDML

I = strcmp(tf, drosTF.names);
tflabel = drosTF.labels(I);
inds = drosTF.chipinds{I};

J = find(sum(droschip.data(:, setdiff(1:15, inds))') == 0);
[foo, I] = sort(-sum(droschip.data(J, :) ~= 0, 2));
genes = droschip.genes(J(I));
