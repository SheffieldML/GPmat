function [states, ll] = viterbiAlign(logTrans, logLike)

% VITERBIALIGN Compute the Viterbi alignment.
% FORMAT
% DESC computes the viterbi alignment given a log transition
% probability matrix and a log emission probability matrix.
% ARG logTrans : the logarithm of the transition probabilities (if
% exponentiated the sum of each row should be one, i.e. the matrix
% gives the transition probabilities from the ith state to the
% other states in the ith row.
% ARG logLike : the logarithm of the emission probabilities. The
% different states are the rows of the matrix.
% RETURN states : the most likely state at each point in the
% sequence.
% RETURN ll : the log likelihood of the most likely state sequence.
%
% SEEALSO :
%
% COPYRIGHT : Neil D. Lawrence, 2006
% 

% MLTOOLS

logProbPlace = zeros(size(logLike));
logProbPlace(:, 1) = logLike(:, 1);
ind(:, 1) = [1:size(logLike, 1)]';
for i = 1:size(logLike, 2)-1
  nextProbs = repmat(logProbPlace(:, i)...
                     +logLike(:, i), 1, size(logTrans, 2)) ...
      + logTrans;
  [logProbPlace(:, i+1), tempInd]= max(nextProbs, [], 2);
  ind(tempInd, i+1) = (1:size(ind, 1))';
end
states = zeros(1, size(logLike, 2));
[ll, states(end)] = max(logProbPlace(:, end));
for i = size(states, 2):-1:2
  states(i-1)= ind(states(i), i-1);
end
