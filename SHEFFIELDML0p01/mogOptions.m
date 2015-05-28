function options = mogOptions(numComp)

% MOGOPTIONS Sets the default options structure for MOG models.
%
%	Description:
%
%	OPTIONS = MOGOPTIONS(NUMCOMPONENTS) sets the default options
%	structure for mixtures of Gaussians models.
%	 Returns:
%	  OPTIONS - structure containing the default options.
%	 Arguments:
%	  NUMCOMPONENTS - number of components in the mixture model.
%	
%
%	See also
%	MOGCREATE


%	Copyright (c) 2006, 2008 Neil D. Lawrence


options.numComponents = numComp;
options.covtype = 'ppca';
% Whether it is an infinite mixture (false by default);
options.isInfinite = false;
options.a1=1;