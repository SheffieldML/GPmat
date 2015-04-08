function model = multiprobitUpdateSites(model, index)

% MULTIPROBITUPDATESITES Update site parameters for multiprobit model.


model = multiprobitUpdateParams(model, index);
model = updateSites(model, index);