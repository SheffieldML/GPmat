function L = noiseLogLikelihood(noise, mu, varsigma, y);

% NOISELOGLIKELIHOOD Return the log-likelihood under the noise model.
%
%	Description:
%
%	NOISELOGLIKELIHOOD(NOISE, MU, VARSIGMA, Y) returns the log
%	likelihood of a data set under the given noise model.
%	 Arguments:
%	  NOISE - the noise structure for which the log likelihood is
%	   required.
%	  MU - input mean locations for the log likelihood.
%	  VARSIGMA - input variance locations for the log likelihood.
%	  Y - target locations for the log likelihood.
%	
%
%	See also
%	NOISEPARAMINIT, NOISELIKELIHOOD


%	Copyright (c) 2004, 2005 Neil D. Lawrence


fhandle = str2func([noise.type 'NoiseLogLikelihood']);
L = fhandle(noise, mu, varsigma, y);


% check if there is a prior over parameters
if isfield(noise, 'priors')
  fhandle = str2func([noise.type 'NoiseExpandParams']);
  params = fhandle(noise);
  for i = 1:length(noise.priors)
    index = noise.priors(i).index;
    L = L + priorLogProb(noise.priors(i), params(index));
  end
end
