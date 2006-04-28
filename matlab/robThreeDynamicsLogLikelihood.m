function ll = robThreeDynamicsLogLikelihood(model)

% ROBTHREEDYNAMICSLOGLIKELIHOOD Give the log likelihood of the robot three dynamics part.
%
% ll = robThreeDynamicsLogLikelihood(model)
%

% Copyright (c) 2006 Neil D. Lawrence
% robThreeDynamicsLogLikelihood.m version 



ll = 0;
for i = 1:size(model.diffX)-1
  covMat = model.lambda*model.diffX(i, :)'*model.diffX(i, :) + ...
      eye(2)*model.sigma2;
  invCovMat = inv(covMat);
  ll = ll -0.5* model.diffX(i+1, :)*invCovMat*model.diffX(i+1, :)';
end

  