function str = optimiDefaultOptimiser

% OPTIMIDEFAULTOPTIMISER Returns the default optimiser to be used.
% FORMAT
% DESC returns the default optimiser, placing the command here
% makes it easier to change the default globally for all toolboxes.
% RETURN str : string which represents the default optimiser
% (currently 'scg', the NETLAB scaled conjugate gradient
% optimiser).
%
% SEEALSO : scg, conjgrad, minimize, optimiDefaultConstraint
%
% COPYRIGHT : Neil D. Lawrence, 2006

% OPTIMI

str = 'scg';

