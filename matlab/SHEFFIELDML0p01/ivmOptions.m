function options = ivmOptions(varargin)

% IVMOPTIONS Return default options for IVM model.
%
%	Description:
%
%	OPTIONS = IVMOPTIONS returns the default options in a structure for
%	a IVM model.
%	 Returns:
%	  OPTIONS - structure containing the default options for the given
%	   approximation type.
%	
%
%	See also
%	IVMCREATE


%	Copyright (c) 2006, 2005 Neil D. Lawrence


% bog-standard kernel.
options.kern = {'rbf', 'bias', 'white'};
options.numActive = 100;
options.noise = 'probit';
options.selectionCriterion = 'entropy';
options.display = 0;
options.kernIters = 100;
options.noiseIters = 0;
options.extIters = 8;

