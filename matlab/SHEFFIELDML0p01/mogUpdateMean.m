function model = mogUpdateMean(model)

% MOGUPDATEMEAN Update the means of an MOG model.
%
%	Description:
%
%	MODEL = MOGUPDATEMEAN(MODEL) updates the mean vectors of a mixtures
%	of Gaussians model.
%	 Returns:
%	  MODEL - the model with updated means.
%	 Arguments:
%	  MODEL - the model which is to be updated.
%	
%
%	See also
%	MOGCREATE, MOGUPDATECOVARIANCE, MOGESTEP


%	Copyright (c) 2006 Neil D. Lawrence


for i = 1:model.m
  if sum(model.posterior(:, i)) ~= 0
    for j = 1:model.d
      model.mean(i, j) = model.posterior(:, i)'*model.Y(:, ...
                                                        j)/sum(model.posterior(:, i));
    end    
  else
    p = exp(model.lnposterior(:, i) - max(model.lnposterior(:, i)));
    for j = 1:model.d
      model.mean(i, j) = p'*model.Y(:, j)/sum(p);
    end
  end
end

