function gpTimeDynamicsDisplay(model, varargin)

% GPTIMEDYNAMICSDISPLAY Display a GP time dynamics model.
%
%	Description:
%
%	GPTIMEDYNAMICSDISPLAY(MODEL, SPACENUM) displays in human readable
%	form the contents of the GP time dynamics model.
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