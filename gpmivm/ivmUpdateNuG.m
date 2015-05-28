function model = ivmUpdateNuG(model, index)

% IVMUPDATENUG Update nu and g parameters associated with noise model.
% FORMAT
% DESC updates the parameter vectors nu and g which depend on the
% noise model used.
% ARG model : the model being updated.
% ARG ind : the index of the point that has been included.
% RETURN model : the model structure with nu and g updated.
% 
% SEEALSO : ivmCreate, ivmAddPoint, ivmEpUpdateM, ivmInit, ivmRemovePoint,
% ivmSelectPoints, 
%
% COPYRIGHT : Neil D. Lawrence, 2004, 2005, 2006

% IVM

if nargin < 2
  index = 1:size(model.y, 1);
end

[model.g(index, :), model.nu(index, :)] = ...
    noiseUpdateNuG(model.noise, ...
                   model.mu(index, :), model.varSigma(index, :), ...
                   model.y(index, :));


if strcmp(model.noise.type, 'cmpnd') & any(model.nu(index, :)< 0) 
  if model.noise.logconcave
    warning('nu less than zero in log concave model.')
  else
    fprintf('nu less than zero\n')
  end
end
