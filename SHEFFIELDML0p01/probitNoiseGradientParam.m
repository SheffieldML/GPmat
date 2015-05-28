function g = probitNoiseGradientParam(noise, mu, varsigma, y)

% PROBITNOISEGRADIENTPARAM Gradient of PROBIT noise's parameters.
%
%	Description:
%
%	G = PROBITNOISEGRADIENTPARAM(NOISE, MU, VARSIGMA, Y) computes the
%	gradient of the log Z of the probit based classification noise model
%	with respect to the of functions with respect to the probit based
%	classification noise's parameters.
%	 Returns:
%	  G - gradients of the log Z with respect to the noise parameters.
%	   The ordering of the vector should match that provided by the
%	   function noiseExtractParam.
%	 Arguments:
%	  NOISE - the noise structure for which the gradients are being
%	   computed.
%	  MU - the input means for which the gradients are being computed.
%	  VARSIGMA - the input variances for which the gradients are being
%	   computed.
%	  Y - the target values for the noise model.
%	
%	
%
%	See also
%	% SEEALSO PROBITNOISEPARAMINIT, PROBITNOISEGRADVALS, NOISEGRADIENTPARAM


%	Copyright (c) 2004, 2005 Neil D. Lawrence



c = y./sqrt(noise.sigma2 + varsigma);
for i = 1:size(mu, 2)
  u(:, i) = c(:, i).*(mu(:, i) + noise.bias(i));
end
g = sum(c.*gradLogCumGaussian(u), 1);