function gpDynamicsDisplay(model, varargin)

% GPDYNAMICSDISPLAY Display a GP dynamics model.
%
%	Description:
%
%	GPDYNAMICSDISPLAY(MODEL, SPACENUM) displays in human readable form
%	the contents of the GP dynamics model.
%	 Arguments:
%	  MODEL - the model structure to be displaced.
%	  SPACENUM - number of spaces to place before displaying model
%	   structure.
%	
%
%	See also
%	GPDISPLAY, GPDYNAMICSCREATE, MODELDISPLAY.


%	Copyright (c) 2006 Neil D. Lawrence


gpDisplay(model, varargin{:});