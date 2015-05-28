function model = leReconstruct(leInfo, y)

% LERECONSTRUCT Reconstruct an LE form component parts.
%
%	Description:
%
%	MODEL = LERECONSTRUCT(LEINFO, Y) takes component parts of an LE
%	model and reconstructs the LE model. The component parts are
%	normally retrieved from a saved file.
%	 Returns:
%	  MODEL - an LE model structure that combines the component parts.
%	 Arguments:
%	  LEINFO - the active set and other information stored in a
%	   structure.
%	  Y - the output target training data for the LE.
%	
%
%	See also
%	LEDECONSTRUCT, LECREATE


%	Copyright (c) 2009 Neil D. Lawrence


  model = leInfo;
  model.Y = y;
    
end
