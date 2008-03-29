function l = priorLogProb(prior, x)

% PRIORLOGPROB Log probability of given prior.
% FORMAT
% DESC wrapper function that computes the log probability of the data under the given prior.
% ARG prior : the prior structure for which the log probability is to be
% computed.
% ARG x : the data for which the log probability is required.
% RETURN l : the log probability of the data.
%
% SEEALSO : priorCreate
%
% COPYRIGHT : Neil D. Lawrence, 2003, 2004
  
% PRIOR

% Compute log prior
fhandle = str2func([prior.type 'PriorLogProb']);
l = fhandle(prior, x);
