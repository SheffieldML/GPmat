function options = kbrOptions(X)

% KBROPTIONS Create a default options structure for the KBR model.
%
%	Description:
%
%	OPTIONS = KBROPTIONS(X) creates a default options structure for the
%	kernel based regression model structure.
%	 Returns:
%	  OPTIONS - the default options structure.
%	 Arguments:
%	  X - the input data for the kernel regression.
%	
%
%	See also
%	KBRCREATE, MODELOPTIONS


%	Copyright (c) 2005, 2006 Neil D. Lawrence


options.kern = 'rbf';
options.X = X;
