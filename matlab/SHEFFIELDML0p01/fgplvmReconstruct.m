function model = fgplvmReconstruct(kern, noise, fgplvmInfo, X, y)

% FGPLVMRECONSTRUCT Reconstruct an FGPLVM from component parts.
%
%	Description:
%
%	MODEL = FGPLVMRECONSTRUCT(KERN, NOISE, FGPLVMINFO, X, Y) takes
%	component parts of an FGPLVM model and reconstructs the FGPLVM
%	model. The component parts are normally retrieved from a saved file.
%	 Returns:
%	  MODEL - an FGPLVM model structure that combines the component
%	   parts.
%	 Arguments:
%	  KERN - a kernel structure for the FGPLVM.
%	  NOISE - a noise structure for the FGPLVM.
%	  FGPLVMINFO - the active set and the inactive set of the FGPLVM as
%	   well as the site parameters, stored in a structure.
%	  X - the input training data for the FGPLVM.
%	  Y - the output target training data for the FGPLVM.
%	
%
%	See also
%	FGPLVMDECONSTRUCT, FGPLVMCREATE, GPRECONSTRUCT


%	Copyright (c) 2009 Neil D. Lawrence


  model = gpReconstruct(kern, noise, fgplvmInfo, X, y);
  model.type = 'fgplvm';
  if isfield(model, 'back') && ~isempty(model.back)
    switch model.back.type
     case 'kbr'
      model.back.X = model.y;
     otherwise
      
    end
  end
  params = fgplvmExtractParam(model);
  model = fgplvmExpandParam(model, params);
end
