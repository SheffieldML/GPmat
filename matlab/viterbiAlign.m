function [states, ll] = viterbiAlign(logTrans, logLike, logPriorProb)

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
% ARG logPriorProb : the prior probability for each start state. If not
% provided it is assumed to be uniform.
% RETURN states : the most likely state at each point in the
% sequence.
% RETURN ll : the log likelihood of the most likely state sequence.
%
% SEEALSO :
%
% COPYRIGHT : Neil D. Lawrence, 2006, 2007
% 

% MLTOOLS

if nargin < 3
  logPriorProb = -log(size(logLike, 1));
end
logProbPlace = zeros(size(logLike));
logProbPlace(:, 1) = logPriorProb + logLike(:, 1);
ind(:, 1) = [1:size(logLike, 1)]';
for i = 1:size(logLike, 2)-1
  % The most likely ways of getting to the next three states.
  nextProbs = repmat(logLike(:, i+1)', size(logTrans, 1), 1) ...
      + repmat(logProbPlace(:, i), 1, size(logTrans, 1)) ...
      + logTrans;
  [lp, tempInd]= max(nextProbs, [], 1);
  logProbPlace(:, i+1) = lp';
  % store which state was most likely to send us here.
  ind(:, i) = tempInd';
end
states = zeros(1, size(logLike, 2));
[ll, states(end)] = max(logProbPlace(:, end));
for i = size(states, 2):-1:2
  states(i-1)= ind(states(i), i-1);
end
