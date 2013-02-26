function str = optimiDefaultOptimiser

% OPTIMIDEFAULTOPTIMISER Returns the default optimiser to be used.
%
%	Description:
%
%	STR = OPTIMIDEFAULTOPTIMISER returns the default optimiser, placing
%	the command here makes it easier to change the default globally for
%	all toolboxes.
%	 Returns:
%	  STR - string which represents the default optimiser (currently
%	   'scg', the NETLAB scaled conjugate gradient optimiser).
%	
%
%	See also
%	SCG, CONJGRAD, MINIMIZE, OPTIMIDEFAULTCONSTRAINT


%	Copyright (c) 2006 Neil D. Lawrence


str = 'scg';

