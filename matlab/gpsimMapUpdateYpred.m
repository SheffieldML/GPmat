function model = gpsimMapUpdateYpred(model)

% GPSIMMAPUPDATEYPRED Update the stored numerical prediction for y.
% FORMAT
% DESC updates the stored numerical prediction for y. 
% ARG model : the model for which the prediction for y is to be
% updated.
% RETURN model : the model with the updated prediction of y.
%
% SEEALSO : gpsimMapCreate, gpsimMapFunctionalGradient, gpsimMapFunctionalUpdateW
%
% COPYRIGHT : Magnus Rattray and Neil D. Lawrence, 2006

% Run the dynamics forward from t=0 for the numerical approximation.
for j = 1:model.numGenes
  % The model prediction up to time zero is simply B/D;
  model.ypred(1:model.times_index(1), j) = model.B(j)/model.D(j);
  integral = 0;
  for i = model.times_index(1)+1:model.numMapPts
    tfs = model.mapt(i);
    % Compute the integral from zero to current time.
    integral = integral+model.g(i)*exp(model.D(j)*tfs)*model.step;
    % Prediction is integral up to current time multiplied by
    % S*exp(-D*time) plus B/D.
    model.ypred(i,j) = model.B(j)/model.D(j) + model.S(j)*exp(-model.D(j)*tfs)*integral;
  end
end
