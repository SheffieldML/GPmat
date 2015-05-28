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
%
% MODIFIED : Pei Gao, 2008  
% SHEFFIELDML

% Run the dynamics forward from t=0 for the numerical approximation.
start = model.times_index(1)+1;
start = 1;
for j = 1:model.numGenes
    % The model prediction up to time zero is simply B/D;
    % model.ypred(1:model.times_index(1), j) = model.B(j)/model.D(j);
    lnintegral = zeros(1, model.numMapPts);
    
    if isfield(model,'ngParam') && model.ngParam
      gInd = j;
    else
      gInd = 1;
    end
    for i = 1:model.numMapPts
      tfs = model.mapt(i);      
      if isfield(model,'isGroupNonlinearity')
        base = model.alpha(j)*exp(-model.D(j)*tfs); 
      else
        base = 0;
      end
      
      % Compute the integral from zero to current time.
      lnintegral(i) = (model.D(j)*tfs)+log(model.step);
        % Prediction is integral up to current time multiplied by
        % S*exp(-D*time) plus B/D.
      model.ypred(i,j) = base + model.B(j)/model.D(j) ...
       + sum(model.g(start:i, gInd)'.*exp(log(model.S(j))-model.D(j)*tfs+lnintegral(1, start:i)));
    end
end
