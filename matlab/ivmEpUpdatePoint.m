function model = ivmEpUpdatePoint(model, i)

% IVMEPUPDATEPOINT Do an EP update of a point.

% IVM

index = find(model.I == i);
if isempty(index)
  error(['Point ' num2str(i) ' is not in active set'])
end

%/~model = ivmUpdateNuG(model, i);
%~/
% Set nu ready for the point removal.
model = ivmDowndateNuG(model, i);
model = ivmEpUpdateM(model, i);
