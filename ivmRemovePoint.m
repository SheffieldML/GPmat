 function model = ivmRemovePoint(model, i)

% IVMREMOVEPOINT Removes a given point from the IVM.
% FORMAT
% DESC removes a given point from the IVM.
% ARG model : the model from which the point is to be removed.
% ARG ind : the index of the point to be removed.
% RETURN model : the model with the point removed.
% 
% SEEALSO : ivmDowndateNuG, ivmDowndateM, ivmDowndateSites,
% ivmSelectPoints
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006

% IVM

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
