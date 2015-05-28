function model = mvuReconstruct(mvuInfo, y)

% MVURECONSTRUCT Reconstruct an MVU form component parts.
%
%	Description:
%
%	MODEL = MVURECONSTRUCT(MVUINFO, Y) takes component parts of an MVU
%	model and reconstructs the MVU model. The component parts are
%	normally retrieved from a saved file.
%	 Returns:
%	  MODEL - an MVU model structure that combines the component parts.
%	 Arguments:
%	  MVUINFO - the active set and other information stored in a
%	   structure.
%	  Y - the output target training data for the MVU.
%	
%
%	See also
%	MVUDECONSTRUCT, MVUCREATE


%	Copyright (c) 2009 Neil D. Lawrence


  model = mvuInfo;
  model.Y = y;
    
end
