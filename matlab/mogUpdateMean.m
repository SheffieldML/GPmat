function model = mogUpdateMean(model)

% MOGUPDATEMEAN Update the means of an MOG model.
% FORMAT
% DESC updates the mean vectors of a mixtures of
% Gaussians model. 
% ARG model : the model which is to be updated.
% RETURN model : the model with updated means.
%
% SEEALSO : mogCreate, mogUpdateCovariance, mogEstep
%
% COPYRIGHT : Neil D. Lawrence, 2006

% MLTOOLS

for i = 1:model.m
  for j = 1:model.d
    model.mean(i, j) = model.posterior(:, i)'*model.Y(:, ...
                                                      j)/sum(model.posterior(:, i));
  end    
end

