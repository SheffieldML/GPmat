function [m, C] = mogMeanCov(model)

% MOGMEANCOV Project a mixture of Gaussians to a low dimensional space.
%
%	Description:
%
%	[M, C] = MOGMEANCOV(MODEL) computes the mean and covariance of the
%	distribution given by a mixture of Gaussians.
%	 Returns:
%	  M - mean of the mixture of Gaussians.
%	  C - covariance of the mixture of Gaussians.
%	 Arguments:
%	  MODEL - model for which mean and covariance is required.
%	
%
%	See also
%	MOGCREATE


%	Copyright (c) 2008 Neil D. Lawrence


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