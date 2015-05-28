function model = ivmUpdateNuG(model, index)

% IVMUPDATENUG Update nu and g parameters associated with noise model.
%
%	Description:
%
%	MODEL = IVMUPDATENUG(MODEL, IND) updates the parameter vectors nu
%	and g which depend on the noise model used.
%	 Returns:
%	  MODEL - the model structure with nu and g updated.
%	 Arguments:
%	  MODEL - the model being updated.
%	  IND - the index of the point that has been included.
%	ivmSelectPoints,
%	
%
%	See also
%	IVMCREATE, IVMADDPOINT, IVMEPUPDATEM, IVMINIT, IVMREMOVEPOINT, 


%	Copyright (c) 2004, 2005, 2006 Neil D. Lawrence


if nargin < 2
  index = 1:size(model.y, 1);
end

[model.g(index, :), model.nu(index, :)] = ...
    noiseUpdateNuG(model.noise, ...
                   model.mu(index, :), model.varSigma(index, :), ...
                   model.y(index, :));


if strcmp(model.noise.type, 'cmpnd') & any(model.nu(index, :)< 0) 
  if model.noise.logconcave
    warning('nu less than zero in log concave model.')
  else
    fprintf('nu less than zero\n')
  end
end
