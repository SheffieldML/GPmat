function [Xout, logTrans, logLike] =  fgplvmViterbiSequence(model, ...
                                                  X, Y, logTrans, logLike);

% FGPLVMVITERBISEQUENCE Viterbi align a latent sequence.
% FORMAT
% DESC finds a set of training points which are the most likely to
% have generated the given data. The Viterbi algorithm is used to
% determine the most likely sequence of training points.
% given an initialisation and an observed sequence in data space.
% ARG model : the model for which the latent sequence is required.
% ARG Y : the observed sequence in data space.
% ARG logTrans : the logarithm of the transition probabilities
% between latent points (if not given, it is computed).
% ARG logLike : the logarithm of the emission probabilities.
% RETURN X : the most likely training latent points to have
% generated the sequence under the HMM approximation.
% RETURN logTrans : the logarithm of the transition probabilities.
% RETURN logLike : the logarithm of the emission probabilities.
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2007
%
% SEEALSO : fgplvmCreate, fgplvmOptimiseSequence, viterbiAlign

% FGPLVM

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
