function model = ivmAddPoint(model, i)

% IVMADDPOINT Add a point into the IVM representation.
% FORMAT
% DESC incorporates the ith point from the data set into the IVM
% model.
% ARG model : the model to which the point is to be added.
% ARG index : the index of the point in the training data which is
% to be added.
% RETURN model : the returned model with the point added in.
%
% SEEALSO : ivmUpdateSites, ivmUpdateM, ivmUpdateNuG, ivmSelectPoint,
% ivmRemovePoint, ivmCreate
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005

% IVM

index = find(model.J == i);
if isempty(index)
  error(['Point ' num2str(i) ' is not in inactive set'])
end

%/~model = ivmUpdateNuG(model, i);
%~/
model = ivmUpdateSites(model, i);
model = ivmUpdateM(model, i);

% Remove point from the non-active set and place in the active.
model.J(index) = [];
model.I = [model.I; i];

model = ivmUpdateNuG(model, model.J);
