function g = fgplvmPointLogLikeGradient(model, x, y)

% FGPLVMPOINTLOGLIKEGRADIENT Log-likelihood gradient for of a point of the GP-LVM.
%
% g = fgplvmPointLogLikeGradient(model, x, y)
%

% Copyright (c) 2006 Neil D. Lawrence
% fgplvmPointLogLikeGradient.m version 1.1



logTwoPi = log(2*pi);
[mu, varSigma] = gpPosteriorMeanVar(model, x);
[dmu, dvs] = gpPosteriorGradMeanVar(model, x);


% For more general models this should be done with the noise
% toolbox (see ivmGradX in the ivm toolbox for more details).
nu = 1./varSigma;
dlnZ_dmu = zeros(size(nu));
for i = 1:model.d
  dlnZ_dmu(:, i) = y(:, i) - mu(:, i);
end
dlnZ_dmu = dlnZ_dmu.*nu;
dlnZ_dvs = 0.5*(dlnZ_dmu.*dlnZ_dmu - nu);

g = dlnZ_dmu*dmu' + dlnZ_dvs*dvs';

if isfield(model, 'prior')
  g = g + priorGradient(model.prior, x);
end
