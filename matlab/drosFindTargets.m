function genes = drosFindTargets(droschip),

% DROSFINDTARGETS returns a sorted list of targets
%
% Usage:
%   genes = drosFindTargets(droschip)
% COPYRIGHT : Antti Honkela, 2008
  
% GPSIM 

[foo, I] = sort(-sum(droschip.data ~= 0, 2));
genes = droschip.genes(I);
