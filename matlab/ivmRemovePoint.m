 function model = ivmRemovePoint(model, i)

% IVMREMOVEPOINT Remove the least informative point from the IVM.


index = find(model.I == i);
if isempty(index)
  error(['Point ' num2str(i) ' is not in active set']);
end

model = downdateM(model, i);

model = feval([model.noise.type 'UpdateParams'], model, model.J);