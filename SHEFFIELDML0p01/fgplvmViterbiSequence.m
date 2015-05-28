function [Xout, logTrans, logLike] =  fgplvmViterbiSequence(model, ...
                                                  X, Y, logTrans, logLike);

% FGPLVMVITERBISEQUENCE Viterbi align a latent sequence.
%
%	Description:
%
%	[X, LOGTRANS, LOGLIKE] = FGPLVMVITERBISEQUENCE(MODEL, Y, LOGTRANS,
%	LOGLIKE) finds a set of training points which are the most likely to
%	have generated the given data. The Viterbi algorithm is used to
%	determine the most likely sequence of training points. given an
%	initialisation and an observed sequence in data space.
%	 Returns:
%	  X - the most likely training latent points to have generated the
%	   sequence under the HMM approximation.
%	  LOGTRANS - the logarithm of the transition probabilities.
%	  LOGLIKE - the logarithm of the emission probabilities.
%	 Arguments:
%	  MODEL - the model for which the latent sequence is required.
%	  Y - the observed sequence in data space.
%	  LOGTRANS - the logarithm of the transition probabilities between
%	   latent points (if not given, it is computed).
%	  LOGLIKE - the logarithm of the emission probabilities.
%	
%
%	See also
%	FGPLVMCREATE, FGPLVMOPTIMISESEQUENCE, VITERBIALIGN


%	Copyright (c) 2006, 2007 Neil D. Lawrence


if nargin < 4
  for i = 1:size(X, 1)
    for j = 1:size(X, 1)
      if model.dynamics.diff
        target = X(j, :) - X(i, :);
      else
        target = X(j, :);
      end
      logTrans(i, j) = modelPointLogLikelihood(model.dynamics, ...
                                               X(i, :), ...
                                               target);
      
    end
  end
  logTransMax = max(logTrans, [], 2);
  logTrans = logTrans - repmat(logTransMax, 1, size(logTrans, 2));
  logTrans = logTrans - repmat(log(sum(exp(logTrans), 2)), 1, size(logTrans, 2));
end


if nargin < 5
  for i = 1:size(Y, 1)
    for j = 1:size(X, 1);
      logLike(j, i) = fgplvmPointLogLikelihood(model, ...
                                          X(j, :), ...
                                          Y(i, :));
    end
  end
end


ind = viterbiAlign(logTrans, logLike);
Xout = X(ind, :);