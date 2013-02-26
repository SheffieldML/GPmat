function model = gpReconstruct(kern, noise, gpInfo, X, y)

% GPRECONSTRUCT Reconstruct an GP form component parts.
%
%	Description:
%
%	MODEL = GPRECONSTRUCT(KERN, NOISE, GPINFO, X, Y) takes component
%	parts of an GP model and reconstructs the GP model. The component
%	parts are normally retrieved from a saved file.
%	 Returns:
%	  MODEL - an GP model structure that combines the component parts.
%	 Arguments:
%	  KERN - a kernel structure for the GP.
%	  NOISE - a noise structure for the GP (currently ignored).
%	  GPINFO - the active set and other information stored in a
%	   structure.
%	  X - the input training data for the GP.
%	  Y - the output target training data for the GP.
%	
%
%	See also
%	GPDECONSTRUCT, GPCREATE


%	Copyright (c) 2007, 2009 Neil D. Lawrence


  model = gpInfo;
  model.X = X;
  model.y = y;
  model.kern = kern;
  if ~isempty(noise)
    model.noise = noise;
  end
  model.m = gpComputeM(model);
  
  if isfield(model, 'computeS') && model.computeS 
    model.S = model.m*model.m';
  end
  params = gpExtractParam(model);
  model = gpExpandParam(model, params);
  
end
