function model = ivmAddPoint(model, i)

% IVMADDPOINT Add a point.

% IVM

index = find(model.J == i);
if isempty(index)
  error(['Point ' num2str(i) ' is not in inactive set'])
end

model = feval([model.noise.type 'UpdateSites'], model, i);
model = updateM(model, i);

% Remove point from the non-active set and place in the active.
model.J(index) = [];
model.I = [model.I; i];
  
model = feval([model.noise.type 'UpdateParams'], model, model.J);

