
% IVMREMOVEPOINT Removes a given point from the IVM.
%
%	Description:
%
%	MODEL = IVMREMOVEPOINT(MODEL, IND) removes a given point from the
%	IVM.
%	 Returns:
%	  MODEL - the model with the point removed.
%	 Arguments:
%	  MODEL - the model from which the point is to be removed.
%	  IND - the index of the point to be removed.
%	ivmSelectPoints
%	
%
%	See also
%	IVMDOWNDATENUG, IVMDOWNDATEM, IVMDOWNDATESITES, 


%	Copyright (c) 2005, 2006 Neil D. Lawrence


index = find(model.I == i);
if isempty(index)
  error(['Point ' num2str(i) ' is not in active set']);
end


model = ivmDowndateNuG(model, i);
model = ivmDowndateM(model, i);
model = ivmDowndateSites(model, i);

model.I(index) = [];
model.J = [model.J i];

model = ivmUpdateNuG(model, i);
