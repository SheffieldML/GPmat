function [m, C] = mogMeanCov(model)
 
% MOGMEANCOV Project a mixture of Gaussians to a low dimensional space.
% FORMAT
% DESC computes the mean and covariance of the distribution given by a
% mixture of Gaussians.
% ARG model : model for which mean and covariance is required.
% RETURN m : mean of the mixture of Gaussians.
% RETURN C : covariance of the mixture of Gaussians.
%
% COPYRIGHT : Neil D. Lawrence, 2008
%
% SEEALSO : mogCreate

% MLTOOLS

m = zeros(1, model.d);
modelSecondMoment = zeros(model.d);
for i = 1:model.m
  m = m + model.prior(i)*model.mean(i, :);
end
switch model.covtype
  case 'ppca'
   for i = 1:model.m
     modelSecondMoment = modelSecondMoment ...
         + model.prior(i)*model.mean(i, :)'*model.mean(i, :) ...
         + model.prior(i)*(model.W{i}*model.W{i}' + eye(model.d)* ...
                           model.sigma2(i));
   end
  case 'spherical'
   for i = 1:model.m
     modelSecondMoment = modelSecondMoment ...
         + model.prior(i)*model.mean(i, :)'*model.mean(i, :) ...
         + model.prior(i)*(eye(model.d)*model.sigma2(i));
   end
end
C = modelSecondMoment - (m'*m);
