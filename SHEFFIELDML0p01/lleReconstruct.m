function model = lleReconstruct(lleInfo, y)

% LLERECONSTRUCT Reconstruct an LLE form component parts.
%
%	Description:
%
%	MODEL = LLERECONSTRUCT(LLEINFO, Y) takes component parts of an LLE
%	model and reconstructs the LLE model. The component parts are
%	normally retrieved from a saved file.
%	 Returns:
%	  MODEL - an LLE model structure that combines the component parts.
%	 Arguments:
%	  LLEINFO - the active set and other information stored in a
%	   structure.
%	  Y - the output target training data for the LLE.
%	
%
%	See also
%	LLEDECONSTRUCT, LLECREATE


%	Copyright (c) 2009 Neil D. Lawrence


  model = lleInfo;
  model.Y = y;
    
end
