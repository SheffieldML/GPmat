function model = dnetReconstruct(mapping, dnetInfo, y)

% DNETRECONSTRUCT Reconstruct an DNET form component parts.
%
%	Description:
%
%	MODEL = DNETRECONSTRUCT(KERN, NOISE, DNETINFO, X, Y) takes component
%	parts of an DNET model and reconstructs the DNET model. The
%	component parts are normally retrieved from a saved file.
%	 Returns:
%	  MODEL - an DNET model structure that combines the component parts.
%	 Arguments:
%	  KERN - a kernel structure for the DNET.
%	  NOISE - a noise structure for the DNET (currently ignored).
%	  DNETINFO - the active set and other information stored in a
%	   structure.
%	  X - the input training data for the DNET.
%	  Y - the output target training data for the DNET.
%	
%
%	See also
%	DNETDECONSTRUCT, DNETCREATE


%	Copyright (c) 2009 Neil D. Lawrence


  model = dnetInfo;
  model.mapping = mapping;
  model.y = y;
  
  params = dnetExtractParam(model);
  model = dnetExpandParam(model, params);
  model = dnetEstep(model);
  
end
