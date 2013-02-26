function options = defaultOptions;

% DEFAULTOPTIONS The default options for optimisation.
%
%	Description:
%
%	OPTIONS = DEFAULTOPTIONS returns a default options vector for
%	optimisation.
%	 Returns:
%	  OPTIONS - the default options vector.
%	
%
%	See also
%	SCG, CONJGRAD, QUASINEW


%	Copyright (c) 2005, 2006 Neil D. Lawrence


options = [0,  1e-4, 1e-4, 1e-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1e-8, 0.1, 0];