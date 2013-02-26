function model = ivmEpUpdatePoint(model, i)

% IVMEPUPDATEPOINT Do an EP update of a point.
%
%	Description:
%
%	MODEL = IVMEPUPDATEPOINT(MODEL, I) performs an EP update for a given
%	point in an IVM model. Be careful when using because repeated EP
%	updates can be unstable.
%	 Returns:
%	  MODEL - the returned model with the given point updated.
%	 Arguments:
%	  MODEL - the model for which the EP update is to apply.
%	  I - the index of the point which is being updated.
%	
%
%	See also
%	IVMDOWNDATENUGG, IVMEPUPDATEM


%	Copyright (c) 2005, 2006 Neil D. Lawrence


index = find(model.I == i);
if isempty(index)
  error(['Point ' num2str(i) ' is not in active set'])
end

% Set nu ready for the point removal.
%model = ivmDowndateNuG(model, i);
model = ivmEpUpdateM(model, i);
