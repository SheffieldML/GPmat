 function model = ivmRemovePoint(model, i)

% IVMREMOVEPOINT Remove the least informative point from the IVM.


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
