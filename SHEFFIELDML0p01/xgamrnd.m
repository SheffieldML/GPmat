function r = xgamrnd(a, b),

% XGAMRND Draw a sample from the gamma distribution.
%
%	Description:
%
%	R = XGAMRND(A, B) Generates a random draw from the gamma
%	distribution. The function uses both Matlab's RAND and RANDN
%	generators, so to seed it you must seed both RAND and RANDN.
%	 Returns:
%	  R - a draw of random numbers from the gamma distribution. If a and
%	   b are matrices, r will have the same shape. If either is scalar, r
%	   will have the same shape as the one that is a matrix.
%	 Arguments:
%	  A - shape parameter of the gamma distribution.
%	  B - rate parameter of the gamma distribuion (inverse scale).
%	
%
%	See also
%	% SEEALSO RAND, RANDN


%	Copyright (c) 2010 Antti Honkela


% Based on Marsaglia and Tsang, "A Simple Method for generating gamma
% variables", ACM Transactions on Mathematical Software 26(3):363-372
% (2000).

% Implemented in Mex.
