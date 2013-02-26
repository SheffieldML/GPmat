function l = priorLogProb(prior, x)

% PRIORLOGPROB Log probability of given prior.
%
%	Description:
%
%	L = PRIORLOGPROB(PRIOR, X) wrapper function that computes the log
%	probability of the data under the given prior.
%	 Returns:
%	  L - the log probability of the data.
%	 Arguments:
%	  PRIOR - the prior structure for which the log probability is to be
%	   computed.
%	  X - the data for which the log probability is required.
%	
%
%	See also
%	PRIORCREATE


%	Copyright (c) 2003, 2004 Neil D. Lawrence


% Compute log prior
fhandle = str2func([prior.type 'PriorLogProb']);
l = fhandle(prior, x);
