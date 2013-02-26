function model = ivmAddPoint(model, i)

% IVMADDPOINT Add a point into the IVM representation.
%
%	Description:
%
%	MODEL = IVMADDPOINT(MODEL, INDEX) incorporates the ith point from
%	the data set into the IVM model.
%	 Returns:
%	  MODEL - the returned model with the point added in.
%	 Arguments:
%	  MODEL - the model to which the point is to be added.
%	  INDEX - the index of the point in the training data which is to be
%	   added.
%	ivmRemovePoint, ivmCreate
%	
%
%	See also
%	IVMUPDATESITES, IVMUPDATEM, IVMUPDATENUG, IVMSELECTPOINT, 


%	Copyright (c) 2004, 2005 Neil D. Lawrence


index = find(model.J == i);
if isempty(index)
  error(['Point ' num2str(i) ' is not in inactive set'])
end

model = ivmUpdateSites(model, i);
model = ivmUpdateM(model, i);

% Remove point from the non-active set and place in the active.
model.J(index) = [];
model.I = [model.I; i];

model = ivmUpdateNuG(model, model.J);
