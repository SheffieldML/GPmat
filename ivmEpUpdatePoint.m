function model = ivmEpUpdatePoint(model, i)

% IVMEPUPDATEPOINT Do an EP update of a point.
% FORMAT
% DESC performs an EP update for a given point in an IVM model. Be
% careful when using because repeated EP updates can be unstable.
% ARG model : the model for which the EP update is to apply.
% ARG i : the index of the point which is being updated.
% RETURN model : the returned model with the given point updated.
%
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
%
% SEEALSO : ivmDowndateNugG, ivmEpUpdateM

% IVM

index = find(model.I == i);
if isempty(index)
  error(['Point ' num2str(i) ' is not in active set'])
end

%/~
model = ivmUpdateNuG(model, i);
%~/
% Set nu ready for the point removal.
%model = ivmDowndateNuG(model, i);
model = ivmEpUpdateM(model, i);
