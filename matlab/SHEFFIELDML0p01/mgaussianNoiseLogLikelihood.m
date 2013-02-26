function L = mgaussianNoiseLogLikelihood(noise, mu, varsigma, y)

% MGAUSSIANNOISELOGLIKELIHOOD Log likelihood of the data under the MGAUSSIAN noise model.
%
%	Description:
%
%	MGAUSSIANNOISELOGLIKELIHOOD(NOISE, MU, VARSIGMA, Y) returns the log
%	likelihood of a data set under the  multiple output Gaussian noise
%	model.
%	 Arguments:
%	  NOISE - the noise structure for which the log likelihood is
%	   required.
%	  MU - input mean locations for the log likelihood.
%	  VARSIGMA - input variance locations for the log likelihood.
%	  Y - target locations for the log likelihood.
%	
%
%	See also
%	MGAUSSIANNOISEPARAMINIT, MGAUSSIANNOISELIKELIHOOD, NOISELOGLIKELIHOOD


%	Copyright (c) 2004, 2005 Neil D. Lawrence



N = size(mu, 1);
D = size(mu, 2);
for i = 1:D
  varsigma(:, i) = varsigma(:, i) + noise.sigma2(i);
  mu(:, i) = mu(:, i) + noise.bias(i);
end
arg = (y - mu);
arg = arg.*arg./varsigma;

% Remove unlabelled data from likelihood.
arg = arg(:);
unlabelled = find(isnan(arg));
arg(unlabelled) = [];
varsigma(unlabelled) = [];
L = - 0.5*sum(sum(log(varsigma))) ...
    - 0.5*sum(sum(arg)) ...
    - 0.5*N*D*log(2*pi);
