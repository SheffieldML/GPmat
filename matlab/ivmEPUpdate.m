function model = ivmEPUpdate(model, i)

% IVMEPUPDATE Do an EP update of a point.


model = ivmRemovePoint(model, i);
model = ivmAddPoint(model, i);
